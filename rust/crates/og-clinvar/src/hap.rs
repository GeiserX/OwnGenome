//! Which haplotype a record's ALT allele sits on, as far as its genotype
//! can say.
//!
//! Two records can be one allele when their signatures are equal: phased
//! (`|`) with the same alleles in the same phase set (FORMAT/PS, or none on
//! both), or unphased with the same GT, `0/1` and `1/0` being one het call.
//! A phased missing allele counts as "no ALT here" (`1|.` matches `1|0`);
//! an unphased one gives no signature, and nor does a GT with no ALT.

use crate::text::{cmp_digits, is_digits, nth, split};

/// A genotype that carries an ALT allele on a haplotype the GT can name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Hap {
    /// The phase set (`.` without FORMAT/PS) and the alleles joined by `|`,
    /// a missing one written `0`.
    Phased { ps: Vec<u8>, alleles: Vec<u8> },
    /// The GT, with a two-allele het written low allele first.
    Unphased(Vec<u8>),
}

/// The value of FORMAT key `key` in the sample column: the last one when
/// the key repeats (`first` false), or the first.
fn format_value<'a>(fmt: &[u8], smp: &'a [u8], key: &[u8], first: bool) -> Option<&'a [u8]> {
    let f = split(fmt, b':');
    let s = split(smp, b':');
    let mut found = None;
    for (i, k) in f.iter().enumerate() {
        if *k == key && i < s.len() {
            found = Some(s[i]);
            if first {
                break;
            }
        }
    }
    found
}

/// The GT of a record (the first GT key), or empty.
pub fn gt_of<'a>(fmt: &[u8], smp: &'a [u8]) -> &'a [u8] {
    format_value(fmt, smp, b"GT", true).unwrap_or(b"")
}

impl Hap {
    /// The signature of a record's FORMAT and sample columns.
    pub fn of(fmt: &[u8], smp: &[u8]) -> Option<Hap> {
        let gt = format_value(fmt, smp, b"GT", false).unwrap_or(b"");
        let ps = format_value(fmt, smp, b"PS", false).unwrap_or(b".");
        if !gt.iter().any(|c| (b'1'..=b'9').contains(c)) {
            return None;
        }
        if gt.contains(&b'|') {
            if gt.contains(&b'/') {
                return None;
            }
            let mut alleles = Vec::with_capacity(gt.len());
            for (i, a) in split(gt, b'|').into_iter().enumerate() {
                if !(is_digits(a) || a == b".") {
                    return None;
                }
                if i > 0 {
                    alleles.push(b'|');
                }
                alleles.extend_from_slice(if a == b"." { b"0" } else { a });
            }
            return Some(Hap::Phased {
                ps: ps.to_vec(),
                alleles,
            });
        }
        if gt.iter().any(|&c| !(c.is_ascii_digit() || c == b'/')) {
            return None;
        }
        let a = split(gt, b'/');
        if a.len() == 2 && cmp_digits(a[0], a[1]).is_gt() {
            return Some(Hap::Unphased([a[1], b"/", a[0]].concat()));
        }
        Some(Hap::Unphased(gt.to_vec()))
    }

    /// The GT of the record joining parts with GTs `g1` and `g2`: the
    /// signature's alleles, with `.` on a haplotype where either part is
    /// missing (not called there, or the other ALT of a split record).
    pub fn joined_gt(&self, g1: &[u8], g2: &[u8]) -> Vec<u8> {
        match self {
            Hap::Unphased(gt) => gt.clone(),
            Hap::Phased { alleles, .. } => {
                let (a, b) = (split(g1, b'|'), split(g2, b'|'));
                let mut out = Vec::with_capacity(alleles.len());
                for (i, m) in split(alleles, b'|').into_iter().enumerate() {
                    if i > 0 {
                        out.push(b'|');
                    }
                    let missing = nth(&a, i + 1) == b"." || nth(&b, i + 1) == b".";
                    out.extend_from_slice(if m == b"0" && missing { b"." } else { m });
                }
                out
            }
        }
    }

    /// The FORMAT and sample columns of a joined record with GT `gt`:
    /// `GT`, or `GT:PS` for a phased pair in a named phase set.
    pub fn format_columns(&self, gt: &[u8]) -> (&'static [u8], Vec<u8>) {
        match self {
            Hap::Phased { ps, .. } if ps.as_slice() != b"." => (b"GT:PS", [gt, b":", ps].concat()),
            _ => (b"GT", gt.to_vec()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sig(fmt: &str, smp: &str) -> Option<Hap> {
        Hap::of(fmt.as_bytes(), smp.as_bytes())
    }

    #[test]
    fn het_calls_written_either_way_are_one_signature() {
        assert_eq!(sig("GT", "1/0"), sig("GT", "0/1"));
        assert_eq!(sig("GT", "1/0"), Some(Hap::Unphased(b"0/1".to_vec())));
        assert_ne!(sig("GT", "0/1"), sig("GT", "1/1"));
        assert_eq!(sig("GT", "1"), Some(Hap::Unphased(b"1".to_vec())));
    }

    #[test]
    fn what_gives_no_signature() {
        for (f, s) in [
            ("GT", "0/0"),
            ("GT", "./1"),
            ("GT", "0|1/1"),
            ("GT", "x|1"),
            ("DP", "1"),
            ("GT", ""),
        ] {
            assert_eq!(sig(f, s), None, "{f} {s}");
        }
    }

    #[test]
    fn phase_sets_and_missing_alleles() {
        assert_eq!(sig("GT:PS", "1|.:7"), sig("GT:PS", "1|0:7"));
        assert_ne!(sig("GT:PS", "1|0:7"), sig("GT:PS", "1|0:8"));
        assert_ne!(sig("GT:PS", "1|0:7"), sig("GT:PS", "0|1:7"));
        let h = sig("GT:PS", "1|.:7").unwrap();
        assert_eq!(h.joined_gt(b"1|.", b"1|0"), b"1|.");
        assert_eq!(h.joined_gt(b"1|0", b"1|0"), b"1|0");
        assert_eq!(h.format_columns(b"1|."), (&b"GT:PS"[..], b"1|.:7".to_vec()));
        let p = sig("GT", "0|1").unwrap();
        assert_eq!(p.format_columns(b"0|1"), (&b"GT"[..], b"0|1".to_vec()));
    }
}
