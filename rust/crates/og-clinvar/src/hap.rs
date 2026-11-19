//! Which haplotype a record's ALT allele sits on, as far as its genotype
//! can say.
//!
//! Two records can be one allele when their signatures are equal: phased
//! (`|`) with the same alleles in the same phase set (FORMAT/PS, or none on
//! both), or unphased with the same GT, `0/1` and `1/0` being one het call.
//! A phased missing allele counts as "no ALT here" (`1|.` matches `1|0`);
//! an unphased one gives no signature, and nor does a GT with no ALT.
//!
//! A combined record is only as sure as its weakest part: [`weakest`] gives
//! the GQ and the DP it carries, so `og clinvar hits` and the reports judge
//! it against the call floor as they judge a single call.

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

/// A call's FORMAT value `key` (GQ or DP) as `og clinvar hits` reads it:
/// the last such key's value when it is a whole number, else none.
pub fn metric<'a>(fmt: &[u8], smp: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    format_value(fmt, smp, key, false).filter(|v| is_digits(v))
}

/// The value of one metric (GQ or DP) on the combined record of parts with
/// values `a` and `b`: the lower of the two; `0` when a part has none while
/// the metric is `measured` (the VCF header defines it) or the other part
/// has it, since nothing measured supports that part; none when the metric
/// is not measured at all, as on every single call of that VCF.
pub fn weakest(a: Option<&[u8]>, b: Option<&[u8]>, measured: bool) -> Option<Vec<u8>> {
    match (a, b) {
        (Some(x), Some(y)) => Some(if cmp_digits(x, y).is_le() { x } else { y }.to_vec()),
        (None, None) if !measured => None,
        _ => Some(b"0".to_vec()),
    }
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
    /// `GT`, then `PS` for a phased pair in a named phase set, then `GQ`
    /// and `DP` when the record has them (from [`weakest`]).
    pub fn format_columns(
        &self,
        gt: &[u8],
        gq: Option<&[u8]>,
        dp: Option<&[u8]>,
    ) -> (Vec<u8>, Vec<u8>) {
        let (mut fmt, mut smp) = (b"GT".to_vec(), gt.to_vec());
        let ps = match self {
            Hap::Phased { ps, .. } if ps.as_slice() != b"." => Some(ps.as_slice()),
            _ => None,
        };
        for (key, value) in [(&b"PS"[..], ps), (b"GQ", gq), (b"DP", dp)] {
            if let Some(v) = value {
                fmt.extend([&b":"[..], key].concat());
                smp.extend([&b":"[..], v].concat());
            }
        }
        (fmt, smp)
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
        assert_eq!(
            h.format_columns(b"1|.", None, None),
            (b"GT:PS".to_vec(), b"1|.:7".to_vec())
        );
        let p = sig("GT", "0|1").unwrap();
        assert_eq!(
            p.format_columns(b"0|1", None, None),
            (b"GT".to_vec(), b"0|1".to_vec())
        );
        assert_eq!(
            h.format_columns(b"1|0", Some(b"3"), Some(b"2")),
            (b"GT:PS:GQ:DP".to_vec(), b"1|0:7:3:2".to_vec())
        );
        assert_eq!(
            p.format_columns(b"0|1", None, Some(b"0")),
            (b"GT:DP".to_vec(), b"0|1:0".to_vec())
        );
    }

    #[test]
    fn a_combined_record_takes_the_weaker_part() {
        let w = |a: Option<&[u8]>, b: Option<&[u8]>, m| weakest(a, b, m);
        assert_eq!(w(Some(b"40"), Some(b"9"), true), Some(b"9".to_vec()));
        assert_eq!(w(Some(b"9"), Some(b"40"), false), Some(b"9".to_vec()));
        assert_eq!(w(Some(b"010"), Some(b"9"), true), Some(b"9".to_vec()));
        assert_eq!(w(Some(b"40"), None, false), Some(b"0".to_vec()));
        assert_eq!(w(None, Some(b"40"), false), Some(b"0".to_vec()));
        assert_eq!(w(None, None, true), Some(b"0".to_vec()));
        assert_eq!(w(None, None, false), None);
    }

    #[test]
    fn metrics_are_whole_numbers_or_nothing() {
        let m = |f: &str, s: &str, k: &str| {
            metric(f.as_bytes(), s.as_bytes(), k.as_bytes()).map(<[u8]>::to_vec)
        };
        assert_eq!(m("GT:GQ:DP", "0/1:3:2", "GQ"), Some(b"3".to_vec()));
        assert_eq!(m("GT:GQ:DP", "0/1:3:2", "DP"), Some(b"2".to_vec()));
        assert_eq!(m("GT:GQ", "0/1:.", "GQ"), None);
        assert_eq!(m("GT:GQ", "0/1:4.5", "GQ"), None);
        assert_eq!(m("GT:GQ", "0/1", "GQ"), None);
        assert_eq!(m("GT", "0/1", "DP"), None);
    }
}
