//! What Python does with a JSON value once `json.load` has read it:
//! `bool(v)`, `str(v)`, `repr(v)`, `v or w` and `d.get(k) or default`, the
//! idioms the summary and the report read their data with.

use std::fmt::{self, Write};

use crate::fmt::repr as float_repr;
use crate::json::Value;

/// Why `str()`, `repr()` or `.get()` raised, or would need what this port
/// does not carry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueError {
    /// AttributeError: `.get()` on a value that is not a dict.
    NotADict,
    /// The `repr()` of a string inside a list or dict, when the string holds
    /// a character above U+00FF. Python escapes it or not by the Unicode
    /// table of its own version (`str.isprintable`); this port carries no
    /// such table, so it refuses instead of guessing. A known limit: `str()`
    /// of a string alone is exact for every character.
    UnicodeRepr,
}

impl ValueError {
    /// The Python exception class, or `Unsupported` for the known limit.
    pub fn exception(self) -> &'static str {
        match self {
            ValueError::NotADict => "AttributeError",
            ValueError::UnicodeRepr => "Unsupported",
        }
    }
}

impl fmt::Display for ValueError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            ValueError::NotADict => "object has no attribute 'get'",
            ValueError::UnicodeRepr => "repr() of a character above U+00FF in a list or dict",
        })
    }
}

impl std::error::Error for ValueError {}

/// `bool(v)`: `None`, `False`, `0`, `0.0`, `-0.0`, `""`, `[]` and `{}` are
/// false; everything else, NaN included, is true.
pub fn truthy(v: &Value) -> bool {
    match v {
        Value::Null => false,
        Value::Bool(b) => *b,
        Value::Int(i) => *i != 0,
        Value::Float(x) => *x != 0.0,
        Value::Str(s) => !s.is_empty(),
        Value::Array(a) => !a.is_empty(),
        Value::Object(o) => !o.is_empty(),
    }
}

/// `a or b`: `a` when it is true, else `b`.
pub fn or<'a>(a: &'a Value, b: &'a Value) -> &'a Value {
    if truthy(a) { a } else { b }
}

/// `d.get(key)`: `None` for a missing key, `Some(&Value::Null)` for a key
/// that holds `null`.
pub fn get<'a>(d: &'a Value, key: &str) -> Result<Option<&'a Value>, ValueError> {
    match d {
        Value::Object(o) => Ok(o.get(key)),
        _ => Err(ValueError::NotADict),
    }
}

/// `d.get(key) or default`: the value when the key is there and the value
/// is true, else `default` (a missing key, `null`, `0`, `""`, `[]`, `{}`).
pub fn get_or<'a>(d: &'a Value, key: &str, default: &'a Value) -> Result<&'a Value, ValueError> {
    Ok(get(d, key)?.map_or(default, |v| or(v, default)))
}

/// `str(v)`, which is also `f"{v}"`: a string as it is, a float by `repr`,
/// `True`, `False`, `None`, and a list or dict by `repr`.
pub fn py_str(v: &Value) -> Result<String, ValueError> {
    match v {
        Value::Str(s) => Ok(s.clone()),
        _ => py_repr(v),
    }
}

/// `repr(v)`: `'text'` quoted and escaped, `[1, 'a']`, `{'k': None}`.
pub fn py_repr(v: &Value) -> Result<String, ValueError> {
    let mut out = String::new();
    write_repr(&mut out, v)?;
    Ok(out)
}

fn write_repr(out: &mut String, v: &Value) -> Result<(), ValueError> {
    match v {
        Value::Null => out.push_str("None"),
        Value::Bool(true) => out.push_str("True"),
        Value::Bool(false) => out.push_str("False"),
        Value::Int(i) => {
            let _ = write!(out, "{i}");
        }
        Value::Float(x) => out.push_str(&float_repr(*x)),
        Value::Str(s) => str_repr(out, s)?,
        Value::Array(a) => {
            out.push('[');
            for (i, x) in a.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                write_repr(out, x)?;
            }
            out.push(']');
        }
        Value::Object(o) => {
            out.push('{');
            for (i, (k, x)) in o.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                str_repr(out, k)?;
                out.push_str(": ");
                write_repr(out, x)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

/// `repr(s)` for a string: single quotes unless the text has a `'` and no
/// `"`; `\\`, the quote, `\t`, `\n`, `\r` escaped; other control
/// characters, U+007F and the unprintable Latin-1 characters (U+0080 to
/// U+00A0 and the soft hyphen U+00AD) as `\xNN`.
fn str_repr(out: &mut String, s: &str) -> Result<(), ValueError> {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            '\0'..='\x1f' | '\x7f'..='\u{a0}' | '\u{ad}' => {
                let _ = write!(out, "\\x{:02x}", u32::from(c));
            }
            '\u{100}'.. => return Err(ValueError::UnicodeRepr),
            c => out.push(c),
        }
    }
    out.push(quote);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::loads;

    fn v(text: &str) -> Value {
        loads(text).unwrap()
    }

    #[test]
    fn truthiness_matches_python() {
        for text in ["null", "false", "0", "0.0", "-0.0", "\"\"", "[]", "{}"] {
            assert!(!truthy(&v(text)), "{text} is false in Python");
        }
        for text in [
            "true",
            "1",
            "-1",
            "NaN",
            "1e-300",
            "\" \"",
            "\"0\"",
            "[0]",
            "{\"\": null}",
        ] {
            assert!(truthy(&v(text)), "{text} is true in Python");
        }
    }

    #[test]
    fn get_or_default() {
        let d = v(r#"{"a": {"x": 1}, "z": 0, "n": null, "e": ""}"#);
        let empty = v("{}");
        assert_eq!(get_or(&d, "a", &empty).unwrap(), &v(r#"{"x": 1}"#));
        for key in ["z", "n", "e", "missing"] {
            assert_eq!(get_or(&d, key, &empty).unwrap(), &empty, "{key}");
        }
        assert_eq!(get(&d, "n").unwrap(), Some(&Value::Null));
        assert_eq!(get(&d, "missing").unwrap(), None);
        assert_eq!(get(&v("[1]"), "a"), Err(ValueError::NotADict));
        assert_eq!(get_or(&Value::Null, "a", &empty), Err(ValueError::NotADict));
        assert_eq!(or(&v("0"), &v("2")), &v("2"));
    }

    #[test]
    fn str_and_repr_match_python() {
        let cases = [
            ("null", "None"),
            ("true", "True"),
            ("1.0", "1.0"),
            ("1e16", "1e+16"),
            ("-0.0", "-0.0"),
            ("NaN", "nan"),
            ("-Infinity", "-inf"),
            (
                "123456789012345678901234567890",
                "123456789012345678901234567890",
            ),
            ("\"é\\n€\"", "é\n€"),
            ("[]", "[]"),
            (
                r#"[1, 2.5, "a", null, true, {"k": [], "it's": "\"q\""}]"#,
                r#"[1, 2.5, 'a', None, True, {'k': [], "it's": '"q"'}]"#,
            ),
            (r#"["it's \"x\""]"#, r#"['it\'s "x"']"#),
            (
                r#"["\\ \t\r\n\u0000\u007f\u0080\u00a0\u00a1\u00ad\u00e9\u00ff"]"#,
                r"['\\ \t\r\n\x00\x7f\x80\xa0¡\xadéÿ']",
            ),
        ];
        for (text, want) in cases {
            assert_eq!(py_str(&v(text)).unwrap(), want, "str({text})");
        }
        assert_eq!(py_repr(&v("\"a\"")).unwrap(), "'a'");
        assert_eq!(py_str(&v(r#"["€"]"#)), Err(ValueError::UnicodeRepr));
        assert_eq!(py_str(&v(r#"{"€": 1}"#)), Err(ValueError::UnicodeRepr));
    }
}
