//! Python's `statistics.median`, `fmean` and `variance`, as CPython 3.15
//! computes them.
//!
//! - [`median`] sorts as Python's `sorted` does (stable, so `-0.0` and
//!   `0.0` keep their input order) and takes the middle value, or
//!   `(a + b) / 2` of the two middle values in float arithmetic.
//!   [`median_int`] is the same over Python ints: the odd median stays an
//!   int, the even one is the correctly rounded float `(a + b) / 2`.
//! - [`fmean`] is `math.fsum(data) / len(data)`: two roundings, as in
//!   CPython, not the correctly rounded mean.
//! - [`variance`] is exact: CPython sums the values and their squares as
//!   fractions and rounds once, at the end. Here every finite float is an
//!   integer times a power of two, so both sums are big integers on a fixed
//!   scale, and the result is one correctly rounded division. A naive
//!   two-pass float formula differs from it in the last digits.

use crate::num::{self, FsumError};
use std::cmp::Ordering;
use std::fmt;

/// Why a statistic raised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// `StatisticsError`: `median` or `fmean` of no data.
    NoData,
    /// `StatisticsError`: `variance` of fewer than two values.
    TooFewPoints,
    /// A NaN in `median`'s data. Python returns a value that depends on
    /// where the NaN sits, because `sorted` cannot order it; this port
    /// refuses instead.
    NotANumber,
    /// `math.fsum` raised inside `fmean`.
    Fsum(FsumError),
    /// `OverflowError`: the exact variance is too large for a float.
    Overflow,
}

impl Error {
    /// The Python exception class (`NotANumber` has none: Python would
    /// return a value).
    pub fn exception(self) -> &'static str {
        match self {
            Error::NoData | Error::TooFewPoints => "StatisticsError",
            Error::NotANumber => "ValueError",
            Error::Fsum(e) => e.exception(),
            Error::Overflow => "OverflowError",
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NoData => f.write_str("no data"),
            Error::TooFewPoints => f.write_str("variance requires at least two data points"),
            Error::NotANumber => f.write_str("NaN in the data of a median"),
            Error::Fsum(e) => e.fmt(f),
            Error::Overflow => f.write_str("integer division result too large for a float"),
        }
    }
}

impl std::error::Error for Error {}

impl From<FsumError> for Error {
    fn from(e: FsumError) -> Self {
        Error::Fsum(e)
    }
}

/// `statistics.median(data)` over floats.
pub fn median<I: IntoIterator<Item = f64>>(data: I) -> Result<f64, Error> {
    let mut v: Vec<f64> = data.into_iter().collect();
    if v.iter().any(|x| x.is_nan()) {
        return Err(Error::NotANumber);
    }
    // Stable, and -0.0 equals 0.0, as `sorted` sees them with `<`.
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let n = v.len();
    match n {
        0 => Err(Error::NoData),
        _ if n % 2 == 1 => Ok(v[n / 2]),
        _ => Ok(v[n / 2]),
    }
}

/// What `statistics.median` returns for ints: an int for an odd count, a
/// float for an even one (`31` or `31.5`, never `31.0` for an odd count).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Median {
    /// The middle value.
    Int(i64),
    /// `(a + b) / 2` of the two middle values.
    Float(f64),
}

impl fmt::Display for Median {
    /// `str()` of the result.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Median::Int(v) => write!(f, "{v}"),
            Median::Float(v) => f.write_str(&crate::fmt::repr(*v)),
        }
    }
}

/// `statistics.median(data)` over ints.
pub fn median_int<I: IntoIterator<Item = i64>>(data: I) -> Result<Median, Error> {
    let mut v: Vec<i64> = data.into_iter().collect();
    let n = v.len();
    if n == 0 {
        return Err(Error::NoData);
    }
    let (low, &mut mid, _) = v.select_nth_unstable(n / 2);
    if n % 2 == 1 {
        return Ok(Median::Int(mid));
    }
    let below = low.iter().copied().max().unwrap_or(mid);
    // Python's int true division is correctly rounded; so is the i128 to
    // f64 conversion, and halving it is exact.
    Ok(Median::Float(
        (i128::from(below) + i128::from(mid)) as f64 / 2.0,
    ))
}

/// `statistics.fmean(data)`: `math.fsum(data) / len(data)`.
pub fn fmean<I: IntoIterator<Item = f64>>(data: I) -> Result<f64, Error> {
    let mut n: u64 = 0;
    let total = num::fsum(data.into_iter().inspect(|_| n += 1))?;
    if n == 0 {
        return Err(Error::NoData);
    }
    Ok(total / n as f64)
}

/// A finite float's magnitude as `m * 2^(shift - 1074)`, `m` below 2^53.
fn decompose(x: f64) -> (u64, usize) {
    let bits = x.abs().to_bits();
    let field = (bits >> 52) as usize;
    let frac = bits & ((1 << 52) - 1);
    if field == 0 {
        (frac, 0)
    } else {
        (frac | 1 << 52, field - 1)
    }
}

/// `statistics.variance(data)`, the sample variance, exact as CPython's.
pub fn variance<I: IntoIterator<Item = f64>>(data: I) -> Result<f64, Error> {
    // Every finite value is ±m * 2^(s - 1074), so with X = Σ ±m 2^s and
    // Y = Σ m² 2^(2s): Σx = X 2^-1074 and Σx² = Y 2^-2148, both exact.
    let (mut pos, mut neg, mut squares) = (Big::default(), Big::default(), Big::default());
    let mut n: u64 = 0;
    // CPython's sx_partials[None]: the non-finite values, summed in order.
    let mut special: Option<f64> = None;
    for x in data {
        n += 1;
        if !x.is_finite() {
            special = Some(special.unwrap_or(0.0) + x);
        } else if special.is_none() {
            let (m, s) = decompose(x);
            let m = u128::from(m);
            if x < 0.0 { &mut neg } else { &mut pos }.add_shifted(m, s);
            squares.add_shifted(m * m, 2 * s);
        }
    }
    if n < 2 {
        return Err(Error::TooFewPoints);
    }
    if let Some(s) = special {
        return Ok(s / (n - 1) as f64);
    }
    // variance = (n Σx² - (Σx)²) / (n (n - 1)) = (nY - X²) 2^-2148 / (n (n - 1)).
    let sum = match pos.cmp(&neg) {
        Ordering::Less => neg.minus(&pos),
        _ => pos.minus(&neg),
    };
    let mut top = squares.times(n).minus(&sum.squared());
    // 192 more bits, so the quotient keeps at least 64 significant bits.
    top.0.splice(0..0, [0; 3]);
    let inexact = top.div_rem(n) != 0;
    let inexact = top.div_rem(n - 1) != 0 || inexact;
    round_scaled(&top, inexact, 2148 + 192)
}

/// `q * 2^-scale` correctly rounded (half to even), `inexact` saying a
/// nonzero remainder sits below `q`'s last bit; `Overflow` past the
/// largest float, as Python's int true division raises.
fn round_scaled(q: &Big, inexact: bool, scale: usize) -> Result<f64, Error> {
    let bits = q.bit_len();
    if bits == 0 {
        return Ok(0.0);
    }
    // Bit j of q weighs 2^(j - scale). Keep 53 bits, or fewer where the
    // value is subnormal and the last bit weighs 2^-1074.
    let low = (bits.saturating_sub(53)).max(scale - 1074);
    let exp = low as i64 - scale as i64;
    if exp > 971 {
        return Err(Error::Overflow);
    }
    let mut m = q.bits_from(low);
    let half = low > 0 && q.bit(low - 1);
    let rest = inexact || (low > 1 && q.any_below(low - 1));
    if half && (rest || m & 1 == 1) {
        m += 1;
    }
    let unit = if exp >= -1022 {
        f64::from_bits(((exp + 1023) as u64) << 52)
    } else {
        f64::from_bits(1 << (exp + 1074))
    };
    // m <= 2^53 is exact, and so is its product with a power of two; only
    // the round-up to 2^1024 overflows.
    let v = m as f64 * unit;
    if v.is_infinite() {
        Err(Error::Overflow)
    } else {
        Ok(v)
    }
}

/// An unsigned integer as little-endian 64-bit limbs, with just the
/// operations [`variance`] needs.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct Big(Vec<u64>);

impl Big {
    fn trim(mut self) -> Self {
        while self.0.last() == Some(&0) {
            self.0.pop();
        }
        self
    }

    /// `self += v << shift`.
    fn add_shifted(&mut self, v: u128, shift: usize) {
        let (at, bit) = (shift / 64, shift % 64);
        let wide = v << bit;
        let high = if bit == 0 {
            0
        } else {
            (v >> (128 - bit)) as u64
        };
        let parts = [wide as u64, (wide >> 64) as u64, high];
        if self.0.len() < at + 3 {
            self.0.resize(at + 3, 0);
        }
        let mut carry = false;
        for (i, limb) in self.0[at..].iter_mut().enumerate() {
            let add = parts.get(i).copied().unwrap_or(0);
            if add == 0 && !carry && i >= parts.len() {
                break;
            }
            let (s, c1) = limb.overflowing_add(add);
            let (s, c2) = s.overflowing_add(u64::from(carry));
            *limb = s;
            carry = c1 || c2;
        }
        if carry {
            self.0.push(1);
        }
    }

    /// `self - other`, for `self >= other`.
    fn minus(&self, other: &Big) -> Big {
        let mut out = self.0.clone();
        let mut borrow = false;
        for (i, limb) in out.iter_mut().enumerate() {
            let sub = other.0.get(i).copied().unwrap_or(0);
            if sub == 0 && !borrow && i >= other.0.len() {
                break;
            }
            let (d, b1) = limb.overflowing_sub(sub);
            let (d, b2) = d.overflowing_sub(u64::from(borrow));
            *limb = d;
            borrow = b1 || b2;
        }
        debug_assert!(!borrow, "minus underflowed");
        Big(out).trim()
    }

    fn times(&self, k: u64) -> Big {
        let mut out = Vec::with_capacity(self.0.len() + 1);
        let mut carry = 0u128;
        for &limb in &self.0 {
            let p = u128::from(limb) * u128::from(k) + carry;
            out.push(p as u64);
            carry = p >> 64;
        }
        out.push(carry as u64);
        Big(out).trim()
    }

    fn squared(&self) -> Big {
        let a = &self.0;
        let mut out = vec![0u64; 2 * a.len()];
        for (i, &x) in a.iter().enumerate() {
            let mut carry = 0u128;
            for (j, &y) in a.iter().enumerate() {
                let t = u128::from(x) * u128::from(y) + u128::from(out[i + j]) + carry;
                out[i + j] = t as u64;
                carry = t >> 64;
            }
            out[i + a.len()] = carry as u64;
        }
        Big(out).trim()
    }

    /// `self /= d`, returning the remainder.
    fn div_rem(&mut self, d: u64) -> u64 {
        let mut rem = 0u128;
        for limb in self.0.iter_mut().rev() {
            let cur = (rem << 64) | u128::from(*limb);
            *limb = (cur / u128::from(d)) as u64;
            rem = cur % u128::from(d);
        }
        *self = std::mem::take(self).trim();
        rem as u64
    }

    fn bit_len(&self) -> usize {
        self.0
            .iter()
            .rposition(|&l| l != 0)
            .map_or(0, |i| 64 * i + 64 - self.0[i].leading_zeros() as usize)
    }

    fn bit(&self, i: usize) -> bool {
        self.0.get(i / 64).is_some_and(|l| l >> (i % 64) & 1 == 1)
    }

    /// Whether any bit below `i` is set.
    fn any_below(&self, i: usize) -> bool {
        let (full, bit) = (i / 64, i % 64);
        self.0[..full.min(self.0.len())].iter().any(|&l| l != 0)
            || (bit > 0 && self.0.get(full).is_some_and(|l| l & ((1 << bit) - 1) != 0))
    }

    /// `self >> i`, which the caller knows fits 64 bits.
    fn bits_from(&self, i: usize) -> u64 {
        let (at, bit) = (i / 64, i % 64);
        let lo = self.0.get(at).copied().unwrap_or(0) >> bit;
        let hi = if bit == 0 {
            0
        } else {
            self.0.get(at + 1).copied().unwrap_or(0) << (64 - bit)
        };
        lo | hi
    }
}

impl PartialOrd for Big {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Big {
    fn cmp(&self, other: &Self) -> Ordering {
        let (a, b) = (self.clone().trim(), other.clone().trim());
        a.0.len()
            .cmp(&b.0.len())
            .then_with(|| a.0.iter().rev().cmp(b.0.iter().rev()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_takes_the_middle_or_the_mean_of_two() {
        assert_eq!(median([3.0, 1.0, 2.0]), Ok(2.0));
        assert_eq!(median([4.0, 1.0, 3.0, 2.0]), Ok(2.5));
        assert_eq!(median([1.0, 2.0]), Ok(1.5));
        assert_eq!(median([1.0, 1.0, 2.0, 2.0]), Ok(1.5));
        assert_eq!(median([7.5]), Ok(7.5));
        assert_eq!(median(Vec::new()), Err(Error::NoData));
        assert_eq!(median([1.0, f64::NAN, 2.0]), Err(Error::NotANumber));
        assert_eq!(median([1e308, 1e308]), Ok(f64::INFINITY));
        assert!(median([f64::INFINITY, f64::NEG_INFINITY]).unwrap().is_nan());
    }

    #[test]
    fn median_keeps_the_input_order_of_signed_zeros() {
        // sorted([0.0, -0.0, 1.0]) == [0.0, -0.0, 1.0]: the middle is -0.0.
        let m = median([0.0, -0.0, 1.0]).unwrap();
        assert!(m == 0.0 && m.is_sign_negative());
        let m = median([-0.0, 0.0, 1.0]).unwrap();
        assert!(m == 0.0 && m.is_sign_positive());
    }

    #[test]
    fn median_int_is_an_int_or_a_rounded_half() {
        assert_eq!(median_int([5, 1, 3]), Ok(Median::Int(3)));
        assert_eq!(median_int([1, 3, 5, 7]), Ok(Median::Float(4.0)));
        assert_eq!(median_int([30, 33, 31, 32]).unwrap().to_string(), "31.5");
        assert_eq!(median_int([31]).unwrap().to_string(), "31");
        assert_eq!(median_int([4, 4]).unwrap().to_string(), "4.0");
        assert_eq!(median_int(Vec::new()), Err(Error::NoData));
        // (2^53 + 1 + 2^53 + 2) / 2 = 2^53 + 1.5, rounded once to 2^53 + 2.
        let a = (1i64 << 53) + 1;
        assert_eq!(
            median_int([a, a + 1]),
            Ok(Median::Float(9007199254740994.0))
        );
        assert_eq!(
            median_int([i64::MAX, i64::MAX]),
            Ok(Median::Float(9223372036854775807.0))
        );
    }

    #[test]
    fn fmean_is_fsum_over_n() {
        assert_eq!(fmean([3.5, 4.0, 5.25]), Ok(4.25));
        assert_eq!(fmean([0.1; 3]), Ok(0.10000000000000002));
        assert_eq!(fmean(Vec::new()), Err(Error::NoData));
        assert_eq!(fmean([1e308, 1e308]), Err(Error::Fsum(FsumError::Overflow)));
        assert_eq!(
            fmean([f64::INFINITY, f64::NEG_INFINITY]),
            Err(Error::Fsum(FsumError::InfMinusInf))
        );
    }

    #[test]
    fn variance_is_exact() {
        // CPython's docstring example.
        let data = [2.75, 1.75, 1.25, 0.25, 0.5, 1.25, 3.5];
        assert_eq!(variance(data), Ok(1.3720238095238095));
        assert_eq!(variance([1.0, 2.0, 3.0, 4.0]), Ok(1.6666666666666667));
        assert_eq!(variance([5.0, 5.0]), Ok(0.0));
        // The naive formula loses everything here; the exact one does not.
        let big = 1e16;
        assert_eq!(variance([big + 2.0, big + 4.0, big + 6.0]), Ok(4.0));
        assert_eq!(variance([0.1, 0.2, 0.3]), Ok(0.009999999999999998));
        assert_eq!(variance([-1.0, 1.0]), Ok(2.0));
        // (5e-324)² / 2 is far below the smallest float.
        assert_eq!(variance([5e-324, 0.0]), Ok(0.0));
    }

    #[test]
    fn variance_edges_as_cpython() {
        assert_eq!(variance([1.0]), Err(Error::TooFewPoints));
        assert_eq!(variance(Vec::new()), Err(Error::TooFewPoints));
        assert_eq!(variance([1.0, f64::INFINITY]), Ok(f64::INFINITY));
        assert_eq!(variance([1.0, f64::INFINITY, 2.0]), Ok(f64::INFINITY));
        assert!(
            variance([f64::INFINITY, f64::NEG_INFINITY])
                .unwrap()
                .is_nan()
        );
        assert!(variance([f64::NAN, 1.0]).unwrap().is_nan());
        assert_eq!(variance([f64::INFINITY]), Err(Error::TooFewPoints));
        // (2e308)^2 / 2 does not fit a float: Python raises OverflowError.
        assert_eq!(variance([1e308, -1e308]), Err(Error::Overflow));
        assert_eq!(variance([1e153, -1e153]), Ok(2e306));
    }

    #[test]
    fn big_arithmetic() {
        let mut a = Big::default();
        a.add_shifted(u128::MAX, 63);
        a.add_shifted(1, 63);
        // (2^128 - 1 + 1) << 63 = 2^191.
        assert_eq!(a.bit_len(), 192);
        assert!(a.bit(191) && !a.any_below(191));
        let mut b = Big(vec![7]).squared().times(3);
        assert_eq!(b, Big(vec![147]));
        assert_eq!(b.div_rem(10), 7);
        assert_eq!(b, Big(vec![14]));
        assert_eq!(Big(vec![0, 1]).minus(&Big(vec![1])), Big(vec![u64::MAX]));
        assert!(Big(vec![0, 1]) > Big(vec![u64::MAX, 0]));
    }
}
