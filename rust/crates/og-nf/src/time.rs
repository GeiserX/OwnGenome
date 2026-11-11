//! The clock formats of the run files, read and written as Python's
//! `datetime` reads and writes them.
//!
//! - [`parse_local`]: a trace time, `YYYY-MM-DD HH:MM:SS[.ffffff]`, naive,
//!   in the local zone of the process (`TZ`, else `/etc/localtime`, else
//!   UTC), as `datetime.strptime(v, fmt).timestamp()`.
//! - [`parse_utc`]: the same fields read as UTC, for the run files' own
//!   formats (`2026-01-02 10:00:00 UTC`, `2026-01-02T10:00:00Z`).
//! - [`iso`]: an epoch as `2026-01-02T10:00:00Z`, as
//!   `datetime.fromtimestamp(e, tz=utc).strftime("%Y-%m-%dT%H:%M:%SZ")`.

use jiff::Timestamp;
use jiff::civil::DateTime;
use jiff::tz::TimeZone;
use pycompat::text::is_space;

/// One piece of a `strptime` format.
#[derive(Clone, Copy)]
pub enum Tok {
    Year,
    Month,
    Day,
    Hour,
    Minute,
    Second,
    Micro,
    /// A literal, matched without regard to ASCII case (strptime compiles
    /// its pattern with IGNORECASE).
    Lit(char),
    /// A space in the format: one or more whitespace characters.
    Space,
}

use Tok::*;

/// `%Y-%m-%d %H:%M:%S.%f`
const TRACE_F: &[Tok] = &[
    Year,
    Lit('-'),
    Month,
    Lit('-'),
    Day,
    Space,
    Hour,
    Lit(':'),
    Minute,
    Lit(':'),
    Second,
    Lit('.'),
    Micro,
];
/// `%Y-%m-%d %H:%M:%S`
const TRACE: &[Tok] = TRACE_F.split_at(11).0;
/// `%Y-%m-%d %H:%M:%S UTC` (run-all.sh's `meta finished_utc`)
pub const RUN_ALL_UTC: &[Tok] = &[
    Year,
    Lit('-'),
    Month,
    Lit('-'),
    Day,
    Space,
    Hour,
    Lit(':'),
    Minute,
    Lit(':'),
    Second,
    Space,
    Lit('U'),
    Lit('T'),
    Lit('C'),
];
/// `%Y-%m-%dT%H:%M:%SZ` (what [`iso`] writes)
pub const ISO_Z: &[Tok] = &[
    Year,
    Lit('-'),
    Month,
    Lit('-'),
    Day,
    Lit('T'),
    Hour,
    Lit(':'),
    Minute,
    Lit(':'),
    Second,
    Lit('Z'),
];

/// The fields a format read: year, month, day, hour, minute, second, µs.
#[derive(Default, Clone, Copy)]
struct Fields([i64; 7]);

/// One alternative of a field: its length and a test of its bytes.
type Alt = (usize, fn(&[u8]) -> bool);

/// The alternatives CPython's `_strptime` tries for a field, in its order:
/// each a length and a test of the digits.
fn alternatives(t: Tok) -> &'static [Alt] {
    fn two(lo: u8, hi: u8) -> impl Fn(&[u8]) -> bool {
        move |d| d[0].is_ascii_digit() && d[1].is_ascii_digit() && (lo..=hi).contains(&num(d))
    }
    fn num(d: &[u8]) -> u8 {
        d.iter().fold(0, |a, c| a * 10 + (c - b'0'))
    }
    fn digit(d: &[u8]) -> bool {
        d[0].is_ascii_digit()
    }
    fn nonzero(d: &[u8]) -> bool {
        matches!(d[0], b'1'..=b'9')
    }
    match t {
        // 1[0-2]|0[1-9]|[1-9]
        Month => &[(2, |d| two(10, 12)(d)), (2, |d| two(1, 9)(d)), (1, nonzero)],
        // 3[0-1]|[1-2]\d|0[1-9]|[1-9]| [1-9]
        Day => &[
            (2, |d| two(30, 31)(d)),
            (2, |d| two(10, 29)(d)),
            (2, |d| two(1, 9)(d)),
            (1, nonzero),
            (2, |d| d[0] == b' ' && nonzero(&d[1..])),
        ],
        // 2[0-3]|[0-1]\d|\d
        Hour => &[(2, |d| two(20, 23)(d)), (2, |d| two(0, 19)(d)), (1, digit)],
        // [0-5]\d|\d
        Minute => &[(2, |d| two(0, 59)(d)), (1, digit)],
        // 6[0-1]|[0-5]\d|\d
        Second => &[(2, |d| two(60, 61)(d)), (2, |d| two(0, 59)(d)), (1, digit)],
        // \d\d\d\d
        Year => &[(4, |d| d.iter().all(u8::is_ascii_digit))],
        // [0-9]{1,6}, greedy
        Micro => &[
            (6, |d| d.iter().all(u8::is_ascii_digit)),
            (5, |d| d.iter().all(u8::is_ascii_digit)),
            (4, |d| d.iter().all(u8::is_ascii_digit)),
            (3, |d| d.iter().all(u8::is_ascii_digit)),
            (2, |d| d.iter().all(u8::is_ascii_digit)),
            (1, digit),
        ],
        Lit(_) | Space => &[],
    }
}

/// Matches `toks` against the whole of `s`, backtracking as the regex does.
fn matches(toks: &[Tok], s: &str, f: &mut Fields) -> bool {
    let Some((&t, rest)) = toks.split_first() else {
        return s.is_empty();
    };
    match t {
        Lit(c) => {
            let mut it = s.chars();
            it.next().is_some_and(|x| x.eq_ignore_ascii_case(&c)) && matches(rest, it.as_str(), f)
        }
        Space => {
            // \s+ is greedy, but what follows never starts with a space.
            let n = s.len() - s.trim_start_matches(is_space).len();
            n > 0 && matches(rest, &s[n..], f)
        }
        _ => {
            let b = s.as_bytes();
            for &(len, ok) in alternatives(t) {
                if b.len() < len || !s.is_char_boundary(len) || !ok(&b[..len]) {
                    continue;
                }
                let digits = std::str::from_utf8(&b[..len]).unwrap_or("").trim_start();
                let mut v: i64 = digits.parse().unwrap_or(0);
                let slot = match t {
                    Year => 0,
                    Month => 1,
                    Day => 2,
                    Hour => 3,
                    Minute => 4,
                    Second => 5,
                    _ => {
                        // "%f": the digits are the leading digits of µs.
                        v *= 10_i64.pow(6 - len as u32);
                        6
                    }
                };
                let keep = f.0[slot];
                f.0[slot] = v;
                if matches(rest, &s[len..], f) {
                    return true;
                }
                f.0[slot] = keep;
            }
            false
        }
    }
}

/// The civil time `s` holds in one of the formats, or `None` when it does
/// not match or is no valid date and time (31 February, second 60).
fn civil(s: &str, formats: &[&[Tok]]) -> Option<(DateTime, i64)> {
    for fmt in formats {
        let mut f = Fields::default();
        if matches(fmt, s, &mut f) {
            let [y, mo, d, h, mi, sec, us] = f.0;
            let dt = DateTime::new(
                i16::try_from(y).ok()?,
                i8::try_from(mo).ok()?,
                i8::try_from(d).ok()?,
                i8::try_from(h).ok()?,
                i8::try_from(mi).ok()?,
                i8::try_from(sec).ok()?,
                0,
            )
            .ok()?;
            return Some((dt, us));
        }
    }
    None
}

/// A trace time (`2026-01-02 10:00:00.000`, or without the fraction) in the
/// zone `tz`, as epoch seconds; `None` when it is not one. A time that the
/// zone skips or repeats is read as Python does with `fold=0`: with the
/// offset in force before the change.
pub fn parse_local(s: &str, tz: &TimeZone) -> Option<f64> {
    let (dt, us) = civil(s, &[TRACE_F, TRACE])?;
    let z = tz.to_ambiguous_zoned(dt).compatible().ok()?;
    Some(z.timestamp().as_second() as f64 + us as f64 / 1e6)
}

/// A time in one of `formats` read as UTC, as epoch seconds.
pub fn parse_utc(s: &str, formats: &[&[Tok]]) -> Option<f64> {
    let (dt, us) = civil(s, formats)?;
    let z = dt.to_zoned(TimeZone::UTC).ok()?;
    Some(z.timestamp().as_second() as f64 + us as f64 / 1e6)
}

/// Epoch seconds as `2026-01-02T10:00:00Z`. The fraction is rounded to
/// microseconds first (half to even), so `9.9999996` is second 10, as in
/// Python. `None` for a time outside the years 1 to 9999.
pub fn iso(epoch: f64) -> Option<String> {
    if !epoch.is_finite() {
        return None;
    }
    let mut secs = epoch.trunc();
    let us = ((epoch - secs) * 1e6).round_ties_even();
    if us >= 1e6 {
        secs += 1.0;
    } else if us < 0.0 {
        secs -= 1.0;
    }
    if secs.abs() > 1e12 {
        return None;
    }
    let ts = Timestamp::from_second(secs as i64).ok()?;
    Some(ts.strftime("%Y-%m-%dT%H:%M:%SZ").to_string())
}

/// A message's `utcTime` (`2026-11-06T10:00:00Z`, `...00.123Z`, `+01:00`,
/// or naive in the zone `tz`) as epoch seconds; `None` for anything else.
/// The shapes Nextflow and run-all.sh write; `datetime.fromisoformat`
/// accepts more.
pub fn parse_iso8601(s: &str, tz: &TimeZone) -> Option<f64> {
    if let Ok(ts) = s.parse::<Timestamp>() {
        return Some(ts.as_second() as f64 + f64::from(ts.subsec_microsecond()) / 1e6);
    }
    let dt: DateTime = s.parse().ok()?;
    let z = tz.to_ambiguous_zoned(dt).compatible().ok()?;
    let ts = z.timestamp();
    Some(ts.as_second() as f64 + f64::from(ts.subsec_microsecond()) / 1e6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_times() {
        let utc = TimeZone::UTC;
        assert_eq!(
            parse_local("2026-01-02 10:00:00.000", &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(
            parse_local("2026-01-02 10:00:00", &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(parse_local("2026-1-2 3:4:5.5", &utc), Some(1_767_323_045.5));
        assert_eq!(
            parse_local("2026-01-02  10:00:00", &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(parse_local("2026-02-30 10:00:00", &utc), None);
        assert_eq!(parse_local("2026-01-02 10:00:60", &utc), None);
        assert_eq!(parse_local("2026-01-02 10:00:00.1234567", &utc), None);
        assert_eq!(parse_local(" 2026-01-02 10:00:00", &utc), None);
        assert_eq!(parse_local("-", &utc), None);
        let madrid = TimeZone::fixed(jiff::tz::offset(1));
        assert_eq!(
            parse_local("2026-01-02 10:00:00", &madrid),
            Some(1_767_344_400.0)
        );
    }

    #[test]
    fn utc_formats() {
        assert_eq!(
            parse_utc("2026-01-02 10:00:00 UTC", &[RUN_ALL_UTC]),
            Some(1_767_348_000.0)
        );
        assert_eq!(
            parse_utc("2026-01-02t10:00:00z", &[ISO_Z]),
            Some(1_767_348_000.0)
        );
        assert_eq!(parse_utc("2026-01-02T10:00:00", &[ISO_Z]), None);
    }

    #[test]
    fn iso_rounds_like_python() {
        assert_eq!(
            iso(1_767_348_000.0).as_deref(),
            Some("2026-01-02T10:00:00Z")
        );
        assert_eq!(
            iso(1_767_348_000.9999996).as_deref(),
            Some("2026-01-02T10:00:01Z")
        );
        assert_eq!(
            iso(1_767_348_000.9).as_deref(),
            Some("2026-01-02T10:00:00Z")
        );
        assert_eq!(iso(-0.5).as_deref(), Some("1969-12-31T23:59:59Z"));
        assert_eq!(iso(f64::NAN), None);
    }

    #[test]
    fn iso8601() {
        let utc = TimeZone::UTC;
        assert_eq!(
            parse_iso8601("2026-01-02T10:00:00Z", &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(
            parse_iso8601("2026-01-02T10:00:00.5Z", &utc),
            Some(1_767_348_000.5)
        );
        assert_eq!(
            parse_iso8601("2026-01-02T11:00:00+01:00", &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(
            parse_iso8601("2026-01-02T10:00:00", &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(parse_iso8601("soon", &utc), None);
    }
}
