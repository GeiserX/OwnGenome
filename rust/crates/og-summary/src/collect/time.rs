//! The minute stamps of summary.json: `datetime.fromtimestamp(ts,
//! tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")`.

/// The UTC minute of a Unix time, as the summary writes it. Python first
/// rounds the time to the microsecond, half to even, so a time a fraction
/// of a microsecond before a minute reads as that minute. `None` for a time
/// Python cannot turn into a date (not finite, or outside years 1 to 9999).
pub fn utc_minute(ts: f64) -> Option<String> {
    if !ts.is_finite() {
        return None;
    }
    let whole = ts.trunc();
    let mut secs = whole as i64;
    let us = ((ts - whole) * 1e6).round_ties_even();
    if us >= 1e6 {
        secs += 1;
    } else if us < 0.0 {
        secs -= 1;
    }
    let days = secs.div_euclid(86_400);
    let in_day = secs.rem_euclid(86_400);
    let (y, m, d) = civil(days);
    if !(1..=9999).contains(&y) {
        return None;
    }
    Some(format!(
        "{y:04}-{m:02}-{d:02} {:02}:{:02} UTC",
        in_day / 3600,
        in_day % 3600 / 60
    ))
}

/// The proleptic Gregorian date of a day count from 1970-01-01 (Howard
/// Hinnant's `civil_from_days`).
fn civil(days: i64) -> (i64, i64, i64) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    (y, m, d)
}

#[cfg(test)]
mod tests {
    use super::utc_minute;

    #[test]
    fn minutes_as_python_writes_them() {
        assert_eq!(utc_minute(0.0).unwrap(), "1970-01-01 00:00 UTC");
        assert_eq!(utc_minute(1_577_836_800.0).unwrap(), "2020-01-01 00:00 UTC");
        assert_eq!(utc_minute(1_709_210_096.5).unwrap(), "2024-02-29 12:34 UTC");
        assert_eq!(utc_minute(-1.0).unwrap(), "1969-12-31 23:59 UTC");
        // Rounded to the microsecond first: half a microsecond before a
        // minute is the minute (half to even: 59.9999995 s rounds up).
        assert_eq!(utc_minute(59.999_999_7).unwrap(), "1970-01-01 00:01 UTC");
        assert_eq!(utc_minute(59.999_99).unwrap(), "1970-01-01 00:00 UTC");
        assert_eq!(utc_minute(f64::NAN), None);
        assert_eq!(utc_minute(f64::INFINITY), None);
        assert_eq!(utc_minute(1e15), None);
    }
}
