//! The read-support limit: the smallest VAF at which a clone shows at least
//! [`MIN_ALT_READS`] reads with the change with probability [`DETECT_POWER`]
//! (binomial). It is an upper bound on what the caller finds: Mutect2 and
//! its filters only take reads away.
//!
//! [`limit`] works at one depth. [`limit_over`] works over a gene's coding
//! bases, each at its own depth: the power at a VAF is the mean over the
//! bases of the power at each base, so a base the run did not read counts
//! as a base where no clone can be seen. At one depth for every base it is
//! [`limit`].
//!
//! The arithmetic is the Python's, step for step, so the numbers agree to
//! the last bit in practice: `C(n, i)` exact (in `u128`) and rounded once to
//! a float, `f ** i` and `(1 - f) ** (n - i)` through `pow`, the terms added
//! with Python's own `sum`, and the same 60-step bisection.

use crate::{DETECT_POWER, Error, MIN_ALT_READS};

/// `C(n, k)`, exact for the small `k` used here (saturating far beyond).
fn comb(n: u64, k: u32) -> u128 {
    let (n, k) = (u128::from(n), u128::from(k));
    if k > n {
        return 0;
    }
    (0..k).fold(1u128, |c, i| c.saturating_mul(n - i) / (i + 1))
}

/// `P(X >= k)` for `X ~ Binomial(n, f)`.
pub fn p_at_least(k: u32, n: u64, f: f64) -> f64 {
    let terms = (0..k).map(|i| {
        let rest = n.saturating_sub(u64::from(i));
        // u128 to f64 rounds to nearest, as Python's int to float does.
        comb(n, i) as f64 * f.powf(f64::from(i)) * (1.0 - f).powf(rest as f64)
    });
    1.0 - pycompat::num::sum(terms)
}

/// The smallest VAF found with probability `power` at `depth` reads, with
/// `k` reads needed; `None` below `k` reads. `depth` is truncated to a whole
/// number of reads, as `int()` does; a depth that is not finite is an error.
pub fn limit(depth: f64, k: u32, power: f64) -> Result<Option<f64>, Error> {
    finite_mean(depth)?;
    let t = depth.trunc();
    if t < f64::from(k) {
        return Ok(None);
    }
    // Depths are far below 2^64; a larger one saturates.
    let n = t as u64;
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if p_at_least(k, n, mid) >= power {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Ok(Some(hi))
}

/// An error for a mean depth that is not finite (a depth file with `nan`
/// or `inf`).
pub fn finite_mean(depth: f64) -> Result<(), Error> {
    if depth.is_finite() {
        Ok(())
    } else {
        Err(Error(format!(
            "a mean depth of {} has no detection limit",
            pycompat::fmt::repr(depth)
        )))
    }
}

/// [`limit`] with the step's own `k` and power.
pub fn detection_limit(depth: f64) -> Result<Option<f64>, Error> {
    limit(depth, MIN_ALT_READS, DETECT_POWER)
}

/// The mean over bases of `P(X >= k)` at VAF `f`. `depths` holds
/// `(depth, bases at that depth)`; a base under `k` reads adds 0.
pub fn power_over(depths: &[(u64, u64)], k: u32, f: f64) -> f64 {
    let (mut hit, mut all) = (0.0f64, 0u64);
    for &(d, n) in depths {
        all += n;
        if d >= u64::from(k) {
            hit += n as f64 * p_at_least(k, d, f);
        }
    }
    if all == 0 { 0.0 } else { hit / all as f64 }
}

/// The smallest VAF at which [`power_over`] reaches `power`, by the
/// same 60-step bisection as [`limit`]; `None` when no VAF under 1 reaches
/// it. The most any VAF reaches is the share of bases with `k` reads or
/// more, so a gene with `1 - power` of its bases or more under `k` reads
/// has no limit.
pub fn limit_over(depths: &[(u64, u64)], k: u32, power: f64) -> Option<f64> {
    let all: u64 = depths.iter().map(|&(_, n)| n).sum();
    let deep: u64 = depths
        .iter()
        .filter(|&&(d, _)| d >= u64::from(k))
        .map(|&(_, n)| n)
        .sum();
    // Only a VAF of 1 reaches a power equal to this share, and near 1 the
    // power of a deep base rounds to exactly 1: test the share itself.
    if all == 0 || deep as f64 / all as f64 <= power {
        return None;
    }
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if power_over(depths, k, mid) >= power {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

/// [`limit_over`] with the step's own `k` and power.
pub fn detection_limit_over(depths: &[(u64, u64)]) -> Option<f64> {
    limit_over(depths, MIN_ALT_READS, DETECT_POWER)
}

/// Whether a clone at VAF `f` shows [`MIN_ALT_READS`] reads with the change
/// at `depth` with probability [`DETECT_POWER`] or more. A base at exactly
/// the depth a limit came from stays readable at that limit: the slack
/// covers the last bit the mean over bases can round away.
pub fn readable_at(depth: u64, f: f64) -> bool {
    depth >= u64::from(MIN_ALT_READS) && p_at_least(MIN_ALT_READS, depth, f) >= DETECT_POWER - 1e-12
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_binomial_coefficients() {
        assert_eq!(comb(30, 0), 1);
        assert_eq!(comb(30, 2), 435);
        assert_eq!(comb(2, 3), 0);
        assert_eq!(comb(6_000_000_000, 2), 17_999_999_997_000_000_000);
    }

    #[test]
    fn limits_python_computes() {
        // detection_limit(d) of bin/chip_report.py (CPython 3.14), repr'd.
        for (d, want) in [
            (30.0, 0.13728764032080973),
            (60.0, 0.06997174201270913),
            (3.0, 0.9283177667225558),
            (45.0, 0.09270293114661994),
            (2000.0, 0.0021382961727698943),
        ] {
            let got = detection_limit(d).unwrap().unwrap();
            assert!((got - want).abs() < 1e-12, "{d}: {got} != {want}");
        }
        assert_eq!(detection_limit(2.9).unwrap(), None);
        assert_eq!(detection_limit(-1.0).unwrap(), None);
        assert!(detection_limit(f64::NAN).is_err());
    }

    #[test]
    fn one_depth_over_every_base_is_the_limit_at_that_depth() {
        for d in [3u64, 15, 30, 33, 60, 2000] {
            let one = detection_limit(d as f64).unwrap().unwrap();
            let over = detection_limit_over(&[(d, 271)]).unwrap();
            assert!((one - over).abs() < 1e-12, "{d}: {one} != {over}");
        }
        assert_eq!(detection_limit_over(&[(2, 100)]), None);
        assert_eq!(detection_limit_over(&[]), None);
    }

    #[test]
    fn an_unread_fifth_leaves_no_limit_at_all() {
        // 80 of 100 bases at 33x: at most 8 in 10 of the bases can show a
        // clone, which only a VAF of 1 reaches. 79 of 100: never.
        assert_eq!(detection_limit_over(&[(33, 80), (0, 20)]), None);
        assert_eq!(detection_limit_over(&[(33, 79), (0, 21)]), None);
        // 90 of 100: a limit, and higher than at 33x everywhere.
        let some = detection_limit_over(&[(33, 90), (0, 10)]).unwrap();
        assert!(some > detection_limit(33.0).unwrap().unwrap(), "{some}");
        assert!((power_over(&[(33, 90), (0, 10)], 3, some) - 0.8).abs() < 1e-6);
    }

    #[test]
    fn readable_at_a_vaf_needs_the_power_at_that_depth() {
        let lim30 = detection_limit(30.0).unwrap().unwrap();
        assert!(readable_at(30, lim30));
        assert!(readable_at(31, lim30));
        assert!(!readable_at(29, lim30));
        assert!(!readable_at(2, 1.0));
    }

    #[test]
    fn red_first_the_limit_needs_k_reads_not_k_plus_one() {
        // P(X >= 3) at the 30x limit is the power; summing i <= k (P(X >= 4))
        // would put it lower.
        let lim = detection_limit(30.0).unwrap().unwrap();
        assert!((p_at_least(3, 30, lim) - 0.8).abs() < 1e-6);
        assert!(p_at_least(4, 30, lim) < 0.7);
    }
}
