//! Python's float and int formatting, byte for byte.
//!
//! Rust's `{:e}` gives the shortest digits that round-trip, which are the
//! digits of Python's `repr`; Rust's `{:.N}` and `{:.Ne}` round the exact
//! binary value half to even, as CPython's dtoa does. What differs is the
//! layout (`1e-05` against `1e-5`, `nan` against `NaN`), and that is what
//! this module writes.

use std::fmt;

/// `nan`, `inf` and `-inf`, the spelling of every Python float format.
fn special(x: f64) -> Option<&'static str> {
    if x.is_nan() {
        Some("nan")
    } else if x == f64::INFINITY {
        Some("inf")
    } else if x == f64::NEG_INFINITY {
        Some("-inf")
    } else {
        None
    }
}

fn sign(x: f64) -> &'static str {
    if x.is_sign_negative() { "-" } else { "" }
}

/// The digits and decimal point of Rust's scientific form of a finite
/// `|x|`: `"1.25e-3"` gives `("125", -2)`, the value `0.125 * 10^-2`.
fn split_sci(s: &str) -> (String, i32) {
    let (mant, exp) = s.split_once('e').expect("scientific form has an exponent");
    let digits = mant.chars().filter(|&c| c != '.').collect();
    let exp: i32 = exp.parse().expect("integer exponent");
    (digits, exp + 1)
}

/// `d[0].d[1..]e±XX`, with at least two exponent digits, as Python writes it.
fn exp_layout(digits: &str, exp: i32) -> String {
    let (first, rest) = digits.split_at(1);
    let dot = if rest.is_empty() { "" } else { "." };
    let es = if exp < 0 { '-' } else { '+' };
    format!("{first}{dot}{rest}e{es}{:02}", exp.unsigned_abs())
}

/// The digits laid out without an exponent, `decpt` digits before the point.
fn fixed_layout(digits: &str, decpt: i32, add_dot_0: bool) -> String {
    let n = digits.len() as i32;
    if decpt <= 0 {
        format!("0.{}{digits}", "0".repeat(decpt.unsigned_abs() as usize))
    } else if decpt >= n {
        let zeros = "0".repeat((decpt - n) as usize);
        let tail = if add_dot_0 { ".0" } else { "" };
        format!("{digits}{zeros}{tail}")
    } else {
        let (int, frac) = digits.split_at(decpt as usize);
        format!("{int}.{frac}")
    }
}

/// The shortest digits that read back as `x` (finite, positive), and the
/// decimal point. Rust and CPython agree on them except when `x` lies
/// exactly halfway between the two nearest candidates: CPython's dtoa then
/// takes the one with the even last digit, Rust may take the other.
fn shortest(x: f64) -> (String, i32) {
    let (digits, decpt) = split_sci(&format!("{x:e}"));
    let Some(exact) = exact_digits(x) else {
        return (digits, decpt);
    };
    let n = digits.len();
    if exact.len() != n + 1 || !exact.ends_with('5') {
        return (digits, decpt);
    }
    // x is the midpoint of `low` and `low + 1` in the last place.
    let low = exact[..n].to_owned();
    let high = increment(&low);
    let parses_back = |d: &str, decpt: i32| format!("0.{d}e{decpt}").parse::<f64>() == Ok(x);
    let (high, high_decpt) = if high.len() > n {
        (high[..n].to_owned(), decpt + 1)
    } else {
        (high, decpt)
    };
    let low_even = (low.as_bytes()[n - 1] - b'0').is_multiple_of(2);
    let mut picks = [(low, decpt), (high, high_decpt)];
    if !low_even {
        picks.swap(0, 1);
    }
    for (d, p) in picks {
        if parses_back(&d, p) {
            let d = d.trim_end_matches('0');
            return (
                if d.is_empty() {
                    "0".to_owned()
                } else {
                    d.to_owned()
                },
                p,
            );
        }
    }
    (digits, decpt)
}

/// A decimal string of digits plus one in the last place (`"19"` -> `"20"`, `"99"` -> `"100"`).
fn increment(d: &str) -> String {
    let mut b = d.as_bytes().to_vec();
    let mut i = b.len();
    loop {
        if i == 0 {
            b.insert(0, b'1');
            break;
        }
        i -= 1;
        if b[i] == b'9' {
            b[i] = b'0';
        } else {
            b[i] += 1;
            break;
        }
    }
    String::from_utf8(b).expect("ascii digits")
}

/// The significant digits of the exact value of a finite positive `x`,
/// when there are at most 18 of them (the only case where a tie between
/// two 17-digit candidates is possible); `None` otherwise.
fn exact_digits(x: f64) -> Option<String> {
    let bits = x.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let frac = bits & ((1 << 52) - 1);
    let (mut m, mut e) = if biased == 0 {
        (frac, -1074)
    } else {
        (frac | (1 << 52), biased - 1075)
    };
    let tz = m.trailing_zeros();
    m >>= tz;
    e += tz as i32;
    let s = if e >= 0 {
        if e < 74 {
            (u128::from(m) << e).to_string()
        } else {
            format!("{x:.0}")
        }
    } else {
        // x = m * 5^k / 10^k; past 5^27 there are more than 18 digits.
        let k = e.unsigned_abs();
        if k > 27 {
            return None;
        }
        (u128::from(m) * 5u128.pow(k)).to_string()
    };
    let s = s.trim_end_matches('0');
    (s.len() <= 18).then(|| s.to_owned())
}

/// `repr(x)`, which is also `str(x)` and `f"{x}"`.
pub fn repr(x: f64) -> String {
    if let Some(s) = special(x) {
        return s.to_owned();
    }
    if x == 0.0 {
        return format!("{}0.0", sign(x));
    }
    let (digits, decpt) = shortest(x.abs());
    let body = if decpt <= -4 || decpt > 16 {
        exp_layout(&digits, decpt - 1)
    } else {
        fixed_layout(&digits, decpt, true)
    };
    format!("{}{body}", sign(x))
}

/// `f"{x:.{prec}f}"`.
pub fn fixed(x: f64, prec: usize) -> String {
    match special(x) {
        Some(s) => s.to_owned(),
        None => format!("{x:.prec$}"),
    }
}

/// `f"{x:.{prec}e}"`.
pub fn sci(x: f64, prec: usize) -> String {
    if let Some(s) = special(x) {
        return s.to_owned();
    }
    let (digits, decpt) = split_sci(&format!("{:.prec$e}", x.abs()));
    let (first, rest) = digits.split_at(1);
    let es = if decpt - 1 < 0 { '-' } else { '+' };
    let dot = if prec == 0 { "" } else { "." };
    format!(
        "{}{first}{dot}{rest}e{es}{:02}",
        sign(x),
        (decpt - 1).unsigned_abs()
    )
}

/// `f"{x:.{prec}g}"`, also `"%.{prec}g" % x`. `f"{x:g}"` is `general(x, 6)`;
/// a precision of 0 counts as 1, as in Python.
pub fn general(x: f64, prec: usize) -> String {
    if let Some(s) = special(x) {
        return s.to_owned();
    }
    if x == 0.0 {
        return format!("{}0", sign(x));
    }
    let p = prec.max(1);
    let (digits, decpt) = split_sci(&format!("{:.*e}", p - 1, x.abs()));
    let digits = digits.trim_end_matches('0');
    let exp = decpt - 1;
    let body = if exp < -4 || exp >= p as i32 {
        exp_layout(digits, exp)
    } else {
        fixed_layout(digits, decpt, false)
    };
    format!("{}{body}", sign(x))
}

/// `f"{x:.{prec}%}"`: `x * 100` in the `f` format, then `%`.
pub fn percent(x: f64, prec: usize) -> String {
    format!("{}%", fixed(x * 100.0, prec))
}

/// Commas every three digits of a run of ASCII digits.
fn group(digits: &str) -> String {
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// `f"{n:,}"` of an int.
pub fn int_grouped(n: i128) -> String {
    let s = if n < 0 { "-" } else { "" };
    format!("{s}{}", group(&n.unsigned_abs().to_string()))
}

/// `f"{x:,.{prec}f}"`.
pub fn fixed_grouped(x: f64, prec: usize) -> String {
    let s = fixed(x, prec);
    let (neg, body) = match s.strip_prefix('-') {
        Some(b) => ("-", b),
        None => ("", s.as_str()),
    };
    if !body.starts_with(|c: char| c.is_ascii_digit()) {
        return s;
    }
    let (int, frac) = match body.split_once('.') {
        Some((i, f)) => (i, Some(f)),
        None => (body, None),
    };
    match frac {
        Some(f) => format!("{neg}{}.{f}", group(int)),
        None => format!("{neg}{}", group(int)),
    }
}

/// `f"{n:0{width}d}"`: the sign first, then zeros up to `width` characters.
pub fn zero_pad(n: i128, width: usize) -> String {
    let digits = n.unsigned_abs().to_string();
    let s = if n < 0 { "-" } else { "" };
    let pad = width.saturating_sub(digits.len() + s.len());
    format!("{s}{}{digits}", "0".repeat(pad))
}

/// `"%-{width}s" % s` and `s.ljust(width)`: width in code points.
pub fn ljust(s: &str, width: usize) -> String {
    let pad = width.saturating_sub(s.chars().count());
    format!("{s}{}", " ".repeat(pad))
}

/// `"%{width}s" % s` and `s.rjust(width)`: width in code points.
pub fn rjust(s: &str, width: usize) -> String {
    let pad = width.saturating_sub(s.chars().count());
    format!("{}{s}", " ".repeat(pad))
}

/// Why Python's `round` raised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoundError {
    /// `round(nan)`: ValueError.
    Nan,
    /// `round(inf)`: OverflowError.
    Infinite,
    /// `round(x, n)` whose result would not fit a float: OverflowError.
    TooLarge,
}

impl RoundError {
    /// The Python exception class.
    pub fn exception(self) -> &'static str {
        match self {
            RoundError::Nan => "ValueError",
            RoundError::Infinite | RoundError::TooLarge => "OverflowError",
        }
    }
}

impl fmt::Display for RoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            RoundError::Nan => "cannot convert float NaN to integer",
            RoundError::Infinite => "cannot convert float infinity to integer",
            RoundError::TooLarge => "rounded value too large to represent",
        })
    }
}

impl std::error::Error for RoundError {}

/// `round(x)`: the nearest integer, ties to even. Python returns an int;
/// this returns it as an integral float, which [`int_str`] prints exactly.
pub fn round_int(x: f64) -> Result<f64, RoundError> {
    if x.is_nan() {
        Err(RoundError::Nan)
    } else if x.is_infinite() {
        Err(RoundError::Infinite)
    } else {
        Ok(x.round_ties_even())
    }
}

/// `str(int(x))` of an integral float: every digit, exactly.
pub fn int_str(x: f64) -> String {
    debug_assert!(x.is_finite() && x.fract() == 0.0);
    if x == 0.0 {
        "0".to_owned()
    } else {
        format!("{x:.0}")
    }
}

/// `round(x, n)` for a float: the exact binary value rounded half to even
/// at `n` decimals (`round(0.35, 1) == 0.3`, `round(2.675, 2) == 2.67`),
/// read back as the nearest float. A negative `n` rounds to tens, hundreds.
pub fn round(x: f64, n: i32) -> Result<f64, RoundError> {
    // CPython's NDIGITS_MAX and NDIGITS_MIN (floatobject.c).
    if !x.is_finite() || n > 323 {
        return Ok(x);
    }
    if n < -308 {
        return Ok(0.0 * x);
    }
    let text = if n >= 0 {
        format!("{x:.0$}", n as usize)
    } else {
        round_integer_digits(x, n.unsigned_abs() as usize)
    };
    let r: f64 = text.parse().expect("a decimal string parses");
    if r.is_infinite() {
        Err(RoundError::TooLarge)
    } else {
        Ok(r)
    }
}

/// `x` rounded half to even to a multiple of `10^k`, as a decimal string
/// with its sign. Works on the exact integer digits so it never rounds twice.
fn round_integer_digits(x: f64, k: usize) -> String {
    let a = x.abs();
    let int = a.trunc();
    let frac_nonzero = a != int;
    let mut digits: Vec<u8> = int_str(int).into_bytes();
    if digits.len() < k + 1 {
        let mut padded = vec![b'0'; k + 1 - digits.len()];
        padded.extend_from_slice(&digits);
        digits = padded;
    }
    let cut = digits.len() - k;
    let (kept, dropped) = digits.split_at(cut);
    let mut kept = kept.to_vec();
    let up = match dropped[0].cmp(&b'5') {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Equal => {
            frac_nonzero
                || dropped[1..].iter().any(|&d| d != b'0')
                || (kept.last().copied().unwrap_or(b'0') - b'0') % 2 == 1
        }
    };
    if up {
        let mut i = kept.len();
        loop {
            if i == 0 {
                kept.insert(0, b'1');
                break;
            }
            i -= 1;
            if kept[i] == b'9' {
                kept[i] = b'0';
            } else {
                kept[i] += 1;
                break;
            }
        }
    }
    let mut s = String::from(sign(x));
    s.push_str(std::str::from_utf8(&kept).expect("ascii digits"));
    s.push_str(&"0".repeat(k));
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repr_matches_python() {
        let cases = [
            (0.1, "0.1"),
            (1.0, "1.0"),
            (-0.0, "-0.0"),
            (0.0, "0.0"),
            (1e16, "1e+16"),
            (1e15, "1000000000000000.0"),
            (1.5e16, "1.5e+16"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (1.2345e-7, "1.2345e-07"),
            (123456.789, "123456.789"),
            (5e-324, "5e-324"),
            (1.7976931348623157e308, "1.7976931348623157e+308"),
            (f64::NAN, "nan"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (2.0f64.powi(60), "1.152921504606847e+18"),
            // Exactly halfway between two 16-digit candidates: the even one.
            (94614632987514.0 + 0.125, "94614632987514.12"),
            (634342696773016.0 + 0.25, "634342696773016.2"),
        ];
        for (x, want) in cases {
            assert_eq!(repr(x), want, "repr({x:?})");
        }
    }

    #[test]
    fn fixed_rounds_the_exact_value_half_even() {
        assert_eq!(fixed(0.125, 2), "0.12");
        assert_eq!(fixed(0.375, 2), "0.38");
        assert_eq!(fixed(2.675, 2), "2.67");
        assert_eq!(fixed(2.5, 0), "2");
        assert_eq!(fixed(-0.04, 1), "-0.0");
        assert_eq!(fixed(f64::NAN, 3), "nan");
        assert_eq!(fixed(-1.0 / 0.0, 3), "-inf");
    }

    #[test]
    fn sci_and_general() {
        assert_eq!(sci(123456.0, 3), "1.235e+05");
        assert_eq!(sci(0.0, 2), "0.00e+00");
        assert_eq!(sci(1.5, 0), "2e+00");
        assert_eq!(general(1234567.0, 6), "1.23457e+06");
        assert_eq!(general(123456.0, 6), "123456");
        assert_eq!(general(0.0001, 6), "0.0001");
        assert_eq!(general(0.00001, 6), "1e-05");
        assert_eq!(general(100.0, 6), "100");
        assert_eq!(general(1.5, 0), "2");
        assert_eq!(general(9.9999996, 6), "10");
        assert_eq!(general(-0.0, 6), "-0");
        assert_eq!(general(0.1 + 0.2, 10), "0.3");
        assert_eq!(general(1e16, 17), "10000000000000000");
    }

    #[test]
    fn percent_grouping_padding() {
        assert_eq!(percent(0.12345, 1), "12.3%");
        assert_eq!(percent(f64::NAN, 1), "nan%");
        assert_eq!(int_grouped(-1234567), "-1,234,567");
        assert_eq!(int_grouped(999), "999");
        assert_eq!(fixed_grouped(1234567.891, 2), "1,234,567.89");
        assert_eq!(fixed_grouped(-999.5, 0), "-1,000");
        assert_eq!(zero_pad(-42, 5), "-0042");
        assert_eq!(zero_pad(42, 1), "42");
        assert_eq!(ljust("é", 3), "é  ");
        assert_eq!(rjust("ab", 3), " ab");
    }

    #[test]
    fn round_is_half_even_on_the_binary_value() {
        assert_eq!(round_int(2.5), Ok(2.0));
        assert_eq!(round_int(3.5), Ok(4.0));
        assert_eq!(round_int(-0.5), Ok(-0.0));
        assert_eq!(round_int(f64::NAN), Err(RoundError::Nan));
        assert_eq!(round_int(f64::INFINITY), Err(RoundError::Infinite));
        assert_eq!(int_str(-0.0), "0");
        assert_eq!(int_str(1e22), "10000000000000000000000");
        assert_eq!(round(0.35, 1), Ok(0.3));
        assert_eq!(round(2.675, 2), Ok(2.67));
        assert_eq!(round(0.125, 2), Ok(0.12));
        assert_eq!(round(1250.0, -2), Ok(1200.0));
        assert_eq!(round(1350.0, -2), Ok(1400.0));
        assert_eq!(round(1250.4, -2), Ok(1300.0));
        assert_eq!(round(5.0, -1), Ok(0.0));
        assert_eq!(round(15.0, -1), Ok(20.0));
        assert_eq!(round(501.0, -3), Ok(1000.0));
        assert!(round(-4.0, -1).unwrap().is_sign_negative());
        assert!(round(-0.04, 1).unwrap().is_sign_negative());
        assert_eq!(
            round(1.7976931348623157e308, -308),
            Err(RoundError::TooLarge)
        );
        assert_eq!(round(1.5, 400), Ok(1.5));
        assert!(round(f64::NAN, 2).unwrap().is_nan());
    }
}
