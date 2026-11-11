//! The detection limit: the smallest VAF a clone needs to be found with
//! probability [`DETECT_POWER`] at a depth, taken as
//! at least [`MIN_ALT_READS`] reads with the change
//! out of that depth (binomial).
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
    if !depth.is_finite() {
        return Err(Error(format!(
            "a mean depth of {} has no detection limit",
            pycompat::fmt::repr(depth)
        )));
    }
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

/// [`limit`] with the step's own `k` and power.
pub fn detection_limit(depth: f64) -> Result<Option<f64>, Error> {
    limit(depth, MIN_ALT_READS, DETECT_POWER)
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
    fn red_first_the_limit_needs_k_reads_not_k_plus_one() {
        // P(X >= 3) at the 30x limit is the power; summing i <= k (P(X >= 4))
        // would put it lower.
        let lim = detection_limit(30.0).unwrap().unwrap();
        assert!((p_at_least(3, 30, lim) - 0.8).abs() < 1e-6);
        assert!(p_at_least(4, 30, lim) < 0.7);
    }
}
