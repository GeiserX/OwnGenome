//! Durations, memory sizes and times as a trace cell or a weblog field
//! holds them: human (`1h 2m 3s`, `18.3 GB`, `2026-01-02 10:00:00.000`) or
//! raw (milliseconds, bytes, epoch milliseconds).

use jiff::tz::TimeZone;
use pycompat::num::parse_float;
use pycompat::text::strip;

/// An empty cell: nothing, or only spaces, or `-`.
pub fn blank(v: Option<&str>) -> bool {
    v.is_none_or(|v| matches!(strip(v, None), "" | "-"))
}

/// Digits with at most one decimal point and digits after it: a raw value.
fn is_number(v: &str) -> bool {
    let (int, frac) = v.split_once('.').unwrap_or((v, "1"));
    !int.is_empty()
        && !frac.is_empty()
        && int.bytes().all(|c| c.is_ascii_digit())
        && frac.bytes().all(|c| c.is_ascii_digit())
}

/// Seconds from a human duration: every `<number><unit>` in it summed
/// (`1h 2m 3s`, `450ms`, `2.5s`, `1d 3h`), units `d`, `h`, `m`, `s`, `ms`.
/// A number glued to another digit or point (`1.2.3s`) counts for nothing.
/// `None` when there is no such pair, or the cell is empty or `-`.
pub fn duration_text(v: &str) -> Option<f64> {
    if v.is_empty() || v == "-" {
        return None;
    }
    let b = v.as_bytes();
    let numeric = |c: u8| c.is_ascii_digit() || c == b'.';
    let (mut total, mut found, mut i) = (0.0, false, 0);
    while i < b.len() {
        if !b[i].is_ascii_digit() || (i > 0 && numeric(b[i - 1])) {
            i += 1;
            continue;
        }
        // The whole run of digits and points from here: it is one number
        // only when it has the shape \d+(\.\d+)?, else nothing in it counts.
        let end = b[i..]
            .iter()
            .position(|&c| !numeric(c))
            .map_or(b.len(), |k| i + k);
        let token = &v[i..end];
        if !is_number(token) {
            i = end;
            continue;
        }
        let rest = v[end..].trim_start_matches(pycompat::text::is_space);
        let (unit, len) = if rest.starts_with("ms") {
            (0.001, 2)
        } else {
            match rest.as_bytes().first() {
                Some(b'd') => (86400.0, 1),
                Some(b'h') => (3600.0, 1),
                Some(b'm') => (60.0, 1),
                Some(b's') => (1.0, 1),
                _ => {
                    i = end;
                    continue;
                }
            }
        };
        total += token.parse::<f64>().unwrap_or(0.0) * unit;
        found = true;
        i = v.len() - rest.len() + len;
    }
    found.then_some(total)
}

/// Seconds from a trace duration: human (see [`duration_text`]) or raw
/// milliseconds.
pub fn duration(v: Option<&str>) -> Option<f64> {
    if blank(v) {
        return None;
    }
    let v = strip(v?, None);
    if is_number(v) {
        return parse_float(v).ok().map(|x| x / 1000.0);
    }
    duration_text(v)
}

/// The bytes of one unit.
fn unit_bytes(u: &str) -> Option<f64> {
    Some(match u.to_ascii_uppercase().as_str() {
        "B" => 1.0,
        "KB" => 1024.0,
        "MB" => 1024f64.powi(2),
        "GB" => 1024f64.powi(3),
        "TB" => 1024f64.powi(4),
        "PB" => 1024f64.powi(5),
        _ => return None,
    })
}

/// Bytes from a memory value: `18.3 GB`, `40 MB`, `64.GB`, or raw bytes,
/// truncated to a whole number as `int(float(...))`. `None` for anything
/// else, a number like `1.2.3` included (the Python this replaces raised
/// there and lost the whole trace).
pub fn memory(v: Option<&str>) -> Option<i128> {
    if blank(v) {
        return None;
    }
    let v = strip(v?, None);
    let x = if is_number(v) {
        parse_float(v).ok()?
    } else {
        let split = v.find(|c: char| !(c.is_ascii_digit() || c == '.'))?;
        let (num, rest) = v.split_at(split);
        let unit = rest.trim_start_matches(pycompat::text::is_space);
        if num.is_empty() || !(1..=2).contains(&unit.len()) || !unit.is_ascii() {
            return None;
        }
        let num = num.parse::<f64>().ok()?;
        num * unit_bytes(unit)?
    };
    x.is_finite().then(|| x.trunc() as i128)
}

/// Epoch seconds from a trace time: `YYYY-MM-DD HH:MM:SS.mmm` in the zone
/// `tz`, or raw epoch milliseconds (0 is no time).
pub fn time(v: Option<&str>, tz: &TimeZone) -> Option<f64> {
    if blank(v) {
        return None;
    }
    let v = strip(v?, None);
    if is_number(v) {
        let ms = parse_float(v).ok()?;
        return (ms > 0.0).then(|| ms / 1000.0);
    }
    crate::time::parse_local(v, tz)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn durations() {
        assert_eq!(duration_text("1h 2m 3s"), Some(3723.0));
        assert_eq!(duration_text("450ms"), Some(0.45));
        assert_eq!(duration_text("2.5s"), Some(2.5));
        assert_eq!(duration_text("1m30s"), Some(90.0));
        assert_eq!(duration_text("1d 1h"), Some(90000.0));
        assert_eq!(duration_text("1.2.3s"), None);
        assert_eq!(duration_text(".5s"), None);
        assert_eq!(duration_text("3 s"), Some(3.0));
        assert_eq!(duration_text("3x 4s"), Some(4.0));
        assert_eq!(duration_text("-"), None);
        assert_eq!(duration_text(""), None);
        assert_eq!(duration(Some("1500")), Some(1.5));
        assert_eq!(duration(Some(" - ")), None);
        assert_eq!(duration(Some("20m 10s")), Some(1210.0));
    }

    #[test]
    fn memories() {
        assert_eq!(
            memory(Some("18.3 GB")),
            Some((18.3 * 1024f64.powi(3)) as i128)
        );
        assert_eq!(memory(Some("40 MB")), Some(40 * 1024 * 1024));
        assert_eq!(memory(Some("64.GB")), Some(64 * 1024 * 1024 * 1024));
        assert_eq!(memory(Some("2147483648")), Some(2_147_483_648));
        assert_eq!(memory(Some("1 kb")), Some(1024));
        assert_eq!(memory(Some("512B")), Some(512));
        assert_eq!(memory(Some("1.2.3 GB")), None);
        assert_eq!(memory(Some(". GB")), None);
        assert_eq!(memory(Some("1 XB")), None);
        assert_eq!(memory(Some("-")), None);
        assert_eq!(memory(None), None);
    }

    #[test]
    fn times() {
        let utc = TimeZone::UTC;
        assert_eq!(time(Some("1767348000000"), &utc), Some(1_767_348_000.0));
        assert_eq!(time(Some("0"), &utc), None);
        assert_eq!(
            time(Some("2026-01-02 10:00:00.000"), &utc),
            Some(1_767_348_000.0)
        );
        assert_eq!(time(Some("soon"), &utc), None);
    }
}
