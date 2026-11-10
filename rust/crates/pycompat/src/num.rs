//! Python's number parsing and float sums.

use std::fmt;

/// `str.isspace()` for one character, the set `str.strip()` removes.
/// Wider than Rust's `char::is_whitespace`: it has `\x1c`-`\x1f`.
pub fn py_isspace(c: char) -> bool {
    matches!(
        c,
        '\t' | '\n'
            | '\x0b'
            | '\x0c'
            | '\r'
            | '\x1c'..='\x1f'
            | ' '
            | '\u{85}'
            | '\u{a0}'
            | '\u{1680}'
            | '\u{2000}'..='\u{200a}'
            | '\u{2028}'
            | '\u{2029}'
            | '\u{202f}'
            | '\u{205f}'
            | '\u{3000}'
    )
}

/// `s.strip()` with no argument.
pub fn py_strip(s: &str) -> &str {
    s.trim_matches(py_isspace)
}

/// The whitespace `float()` and `int()` skip around a number: Unicode
/// whitespace, but of ASCII only C's six (`\x1c`-`\x1f` are refused).
fn number_space(c: char) -> bool {
    py_isspace(c) && !matches!(c, '\x1c'..='\x1f')
}

/// Why `float(s)` or `int(s)` raised.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    /// ValueError: not a number Python reads.
    Invalid,
    /// An int beyond the 128-bit range: Python has no limit, this port has.
    TooLarge,
}

impl ParseError {
    /// The Python exception class (`TooLarge` has none: Python would return the int).
    pub fn exception(&self) -> &'static str {
        match self {
            ParseError::Invalid => "ValueError",
            ParseError::TooLarge => "OverflowError",
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            ParseError::Invalid => "invalid literal",
            ParseError::TooLarge => "integer does not fit 128 bits",
        })
    }
}

impl std::error::Error for ParseError {}

/// Removes underscores that sit between two ASCII digits, as Python allows
/// in a number literal; any other underscore is an error.
fn strip_underscores(s: &str) -> Result<String, ParseError> {
    let b = s.as_bytes();
    let mut out = String::with_capacity(s.len());
    for (i, &c) in b.iter().enumerate() {
        if c == b'_' {
            let before = i > 0 && b[i - 1].is_ascii_digit();
            let after = b.get(i + 1).is_some_and(u8::is_ascii_digit);
            if !(before && after) {
                return Err(ParseError::Invalid);
            }
        } else {
            out.push(c as char);
        }
    }
    Ok(out)
}

/// True when `s` (no sign, no underscores) is a decimal float literal:
/// digits with an optional point, at least one digit, an optional exponent.
fn is_decimal(s: &str) -> bool {
    let b = s.as_bytes();
    let mut i = 0;
    let int = b.iter().take_while(|c| c.is_ascii_digit()).count();
    i += int;
    let mut frac = 0;
    if b.get(i) == Some(&b'.') {
        i += 1;
        frac = b[i..].iter().take_while(|c| c.is_ascii_digit()).count();
        i += frac;
    }
    if int + frac == 0 {
        return false;
    }
    if matches!(b.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(b.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let exp = b[i..].iter().take_while(|c| c.is_ascii_digit()).count();
        if exp == 0 {
            return false;
        }
        i += exp;
    }
    i == b.len()
}

/// `float(s)` for a str: surrounding whitespace, a sign, underscores between
/// digits, `inf`, `infinity` and `nan` in any case; `1e400` is `inf`.
/// Non-ASCII digits, which Python also reads, are refused.
pub fn parse_float(s: &str) -> Result<f64, ParseError> {
    let t = s.trim_matches(number_space);
    if !t.is_ascii() {
        return Err(ParseError::Invalid);
    }
    let (neg, body) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let lower = body.to_ascii_lowercase();
    let v = match lower.as_str() {
        "inf" | "infinity" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let clean = strip_underscores(body)?;
            if !is_decimal(&clean) {
                return Err(ParseError::Invalid);
            }
            clean.parse::<f64>().map_err(|_| ParseError::Invalid)?
        }
    };
    Ok(if neg { -v } else { v })
}

/// `int(s)` for a str in base 10: surrounding whitespace, a sign, leading
/// zeros and underscores between digits. Non-ASCII digits, which Python also
/// reads, are refused, and so is a value beyond `i128`.
pub fn parse_int(s: &str) -> Result<i128, ParseError> {
    let t = s.trim_matches(number_space);
    let (neg, body) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let clean = strip_underscores(body)?;
    if clean.is_empty() || !clean.bytes().all(|c| c.is_ascii_digit()) {
        return Err(ParseError::Invalid);
    }
    let mut v: i128 = 0;
    for d in clean.bytes() {
        let d = i128::from(d - b'0');
        v = v
            .checked_mul(10)
            .and_then(|v| {
                if neg {
                    v.checked_sub(d)
                } else {
                    v.checked_add(d)
                }
            })
            .ok_or(ParseError::TooLarge)?;
    }
    Ok(v)
}

/// Why `math.fsum` raised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FsumError {
    /// OverflowError: finite inputs whose partial sum overflowed.
    Overflow,
    /// ValueError: both `inf` and `-inf` among the inputs.
    InfMinusInf,
}

impl FsumError {
    /// The Python exception class.
    pub fn exception(self) -> &'static str {
        match self {
            FsumError::Overflow => "OverflowError",
            FsumError::InfMinusInf => "ValueError",
        }
    }
}

impl fmt::Display for FsumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            FsumError::Overflow => "intermediate overflow in fsum",
            FsumError::InfMinusInf => "-inf + inf in fsum",
        })
    }
}

impl std::error::Error for FsumError {}

/// `math.fsum(xs)`: the correctly rounded sum (Shewchuk's partials, with
/// CPython's half-even fix-up across partials and its inf/nan rules).
pub fn fsum<I: IntoIterator<Item = f64>>(xs: I) -> Result<f64, FsumError> {
    let mut p: Vec<f64> = Vec::new();
    let mut special = 0.0f64;
    let mut inf_sum = 0.0f64;
    for xsave in xs {
        let mut x = xsave;
        let mut i = 0;
        for j in 0..p.len() {
            let mut y = p[j];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            let lo = y - (hi - x);
            if lo != 0.0 {
                p[i] = lo;
                i += 1;
            }
            x = hi;
        }
        p.truncate(i);
        if x != 0.0 {
            if !x.is_finite() {
                if xsave.is_finite() {
                    return Err(FsumError::Overflow);
                }
                if xsave.is_infinite() {
                    inf_sum += xsave;
                }
                special += xsave;
                p.clear();
            } else {
                p.push(x);
            }
        }
    }
    if special != 0.0 {
        return if inf_sum.is_nan() {
            Err(FsumError::InfMinusInf)
        } else {
            Ok(special)
        };
    }
    let mut n = p.len();
    if n == 0 {
        return Ok(0.0);
    }
    n -= 1;
    let mut hi = p[n];
    let mut lo = 0.0;
    while n > 0 {
        let x = hi;
        n -= 1;
        let y = p[n];
        hi = x + y;
        lo = y - (hi - x);
        if lo != 0.0 {
            break;
        }
    }
    if n > 0 && ((lo < 0.0 && p[n - 1] < 0.0) || (lo > 0.0 && p[n - 1] > 0.0)) {
        let y = lo * 2.0;
        let x = hi + y;
        if y == x - hi {
            hi = x;
        }
    }
    Ok(hi)
}

/// The builtin `sum(xs)` of floats with start 0 (Python 3.12 and later):
/// Neumaier's compensated sum.
pub fn sum<I: IntoIterator<Item = f64>>(xs: I) -> f64 {
    let mut s = 0.0f64;
    let mut c = 0.0f64;
    for x in xs {
        let t = s + x;
        if s.abs() >= x.abs() {
            c += (s - t) + x;
        } else {
            c += (x - t) + s;
        }
        s = t;
    }
    if c != 0.0 && c.is_finite() {
        s += c;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_parsing() {
        assert_eq!(parse_float(" 1_000.5\n"), Ok(1000.5));
        assert_eq!(parse_float("+.5"), Ok(0.5));
        assert_eq!(parse_float("1."), Ok(1.0));
        assert_eq!(parse_float("1e400"), Ok(f64::INFINITY));
        assert_eq!(parse_float("-Infinity"), Ok(f64::NEG_INFINITY));
        assert!(parse_float("nAn").unwrap().is_nan());
        assert_eq!(parse_float("\u{85}2\u{3000}"), Ok(2.0));
        assert_eq!(parse_float("\x1c2"), Err(ParseError::Invalid));
        for bad in [
            "", ".", "1_", "_1", "1__0", "1_.5", "1e", "e5", "1.5.2", "inf1", "0x10", "1 2",
        ] {
            assert_eq!(parse_float(bad), Err(ParseError::Invalid), "{bad:?}");
        }
    }

    #[test]
    fn int_parsing() {
        assert_eq!(parse_int(" -0012 "), Ok(-12));
        assert_eq!(parse_int("1_000"), Ok(1000));
        assert_eq!(parse_int("+7"), Ok(7));
        for bad in ["", "-", "1.0", "1_", "1e3", " "] {
            assert_eq!(parse_int(bad), Err(ParseError::Invalid), "{bad:?}");
        }
        assert_eq!(parse_int(&"9".repeat(40)), Err(ParseError::TooLarge));
        assert_eq!(parse_int(&i128::MIN.to_string()), Ok(i128::MIN));
    }

    #[test]
    fn fsum_is_exact() {
        assert_eq!(fsum([0.1; 10]), Ok(1.0));
        assert_eq!(fsum([1e-16, 1.0, 1e16]), Ok(10000000000000002.0));
        assert_eq!(fsum([1e308, 1e308]), Err(FsumError::Overflow));
        assert_eq!(
            fsum([f64::INFINITY, f64::NEG_INFINITY]),
            Err(FsumError::InfMinusInf)
        );
        assert_eq!(fsum([f64::INFINITY, 1.0]), Ok(f64::INFINITY));
        assert!(fsum([f64::NAN, 1.0]).unwrap().is_nan());
        assert_eq!(fsum(Vec::new()), Ok(0.0));
    }

    #[test]
    fn builtin_sum_is_compensated() {
        assert_eq!(sum([0.1; 10]), 1.0);
        assert_eq!(sum([1e100, 1.0, -1e100, 1.0]), 2.0);
        assert!(!sum([-0.0]).is_sign_negative());
        assert_eq!(sum([f64::INFINITY, 1.0]), f64::INFINITY);
    }
}
