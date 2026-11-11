//! Byte-level field handling with the meaning step 06's tables have always
//! had: tab-separated columns, `|`, `:`, `,` and `;` lists, and the rule
//! that two fields which both read as numbers compare as numbers.

use std::cmp::Ordering;

/// The parts of `s` between `sep` bytes; none for an empty `s`.
pub fn split(s: &[u8], sep: u8) -> Vec<&[u8]> {
    if s.is_empty() {
        Vec::new()
    } else {
        s.split(move |&c| c == sep).collect()
    }
}

/// The 1-based part `i` of `parts`, or empty when there is none.
pub fn nth<'a>(parts: &[&'a [u8]], i: usize) -> &'a [u8] {
    i.checked_sub(1)
        .and_then(|k| parts.get(k))
        .copied()
        .unwrap_or(b"")
}

/// True for one or more ASCII digits and nothing else.
pub fn is_digits(s: &[u8]) -> bool {
    !s.is_empty() && s.iter().all(u8::is_ascii_digit)
}

/// Compares two strings of digits by the numbers they write (leading zeros
/// do not count; an empty string is 0).
pub fn cmp_digits(a: &[u8], b: &[u8]) -> Ordering {
    let strip = |s: &[u8]| -> usize { s.iter().position(|&c| c != b'0').unwrap_or(s.len()) };
    let (a, b) = (&a[strip(a)..], &b[strip(b)..]);
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

/// The number a field reads as when the whole field is one: blanks around
/// it, a sign, digits with an optional fraction, an optional exponent.
pub fn number(s: &[u8]) -> Option<f64> {
    let t = s.trim_ascii();
    let mut i = 0;
    if matches!(t.first(), Some(b'+' | b'-')) {
        i += 1;
    }
    let int = t[i..].iter().take_while(|c| c.is_ascii_digit()).count();
    i += int;
    let mut frac = 0;
    if t.get(i) == Some(&b'.') {
        i += 1;
        frac = t[i..].iter().take_while(|c| c.is_ascii_digit()).count();
        i += frac;
    }
    if int + frac == 0 {
        return None;
    }
    if matches!(t.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(t.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        let exp = t[j..].iter().take_while(|c| c.is_ascii_digit()).count();
        if exp > 0 {
            i = j + exp;
        }
    }
    if i != t.len() {
        return None;
    }
    std::str::from_utf8(t).ok()?.parse().ok()
}

/// Two fields are the same when both read as the same number, or else when
/// their bytes are equal: `1` and `01` name one position.
pub fn same(a: &[u8], b: &[u8]) -> bool {
    match (number(a), number(b)) {
        (Some(x), Some(y)) => x == y,
        _ => a == b,
    }
}

/// A whole number of positions, written as digits only.
pub fn whole(s: &[u8]) -> Option<u64> {
    if !is_digits(s) {
        return None;
    }
    std::str::from_utf8(s).ok()?.parse().ok()
}

/// The value of the first `KEY=value` item of a `;` list, or empty.
pub fn info_get<'a>(info: &'a [u8], key: &[u8]) -> &'a [u8] {
    info.split(|&c| c == b';')
        .find_map(|kv| {
            kv.strip_prefix(key)
                .and_then(|rest| rest.strip_prefix(b"="))
        })
        .unwrap_or(b"")
}

/// `info` with `item` added: the item alone when INFO is empty or `.`.
pub fn info_add(info: &[u8], item: &[u8]) -> Vec<u8> {
    if info.is_empty() || info == b"." {
        item.to_vec()
    } else {
        [info, b";", item].concat()
    }
}

/// `line`'s tab-separated columns with column `k` (1-based) set to `value`;
/// columns up to `k` are added empty when the line has fewer.
pub fn set_column(line: &[u8], k: usize, value: &[u8]) -> Vec<u8> {
    let mut cols = split(line, b'\t');
    while cols.len() < k {
        cols.push(b"");
    }
    cols[k - 1] = value;
    cols.join(&b'\t')
}

/// `s` in capitals (ASCII only).
pub fn upper(s: &[u8]) -> Vec<u8> {
    s.to_ascii_uppercase()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_an_empty_string_gives_no_parts() {
        assert!(split(b"", b'|').is_empty());
        assert_eq!(split(b"a||b", b'|'), vec![&b"a"[..], b"", b"b"]);
        assert_eq!(nth(&split(b"a:b", b':'), 3), b"");
    }

    #[test]
    fn numbers_compare_by_value_and_words_by_bytes() {
        assert!(same(b"1", b"01"));
        assert!(same(b"100", b"1e2"));
        assert!(!same(b"chr1", b"chr01"));
        assert!(!same(b"1", b"1a"));
        assert_eq!(number(b" 2.5 "), Some(2.5));
        assert_eq!(number(b"."), None);
        assert_eq!(number(b"1e"), None);
        assert_eq!(cmp_digits(b"010", b"9"), Ordering::Greater);
        assert_eq!(cmp_digits(b"", b"0"), Ordering::Equal);
    }

    #[test]
    fn info_items_are_read_and_added() {
        assert_eq!(info_get(b"A=1;PGP_X=2;PGP_X=3", b"PGP_X"), b"2");
        assert_eq!(info_get(b"PGP_XY=1", b"PGP_X"), b"");
        assert_eq!(info_add(b".", b"T=1"), b"T=1");
        assert_eq!(info_add(b"A=1", b"T=1"), b"A=1;T=1");
        assert_eq!(set_column(b"a\tb", 4, b"x"), b"a\tb\t\tx");
    }
}
