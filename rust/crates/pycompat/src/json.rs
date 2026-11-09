//! Python's `json` module: `dumps` byte for byte, and `loads` with its
//! leniency (`NaN`, `Infinity`, the last of a duplicate key wins and keeps
//! the first one's place).

use indexmap::IndexMap;
use std::fmt::{self, Write};

use crate::fmt::repr;

/// A JSON value as Python holds it: ints and floats apart, objects in
/// insertion order.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    /// A Python int. Python has no limit; this port stops at 128 bits.
    Int(i128),
    Float(f64),
    Str(String),
    Array(Vec<Value>),
    Object(IndexMap<String, Value>),
}

/// The `json.dumps` keyword arguments this port supports.
#[derive(Debug, Clone)]
pub struct Dump<'a> {
    /// `indent=N`; `None` writes one line.
    pub indent: Option<usize>,
    /// `separators=(item, key)`; `None` is Python's default, which depends
    /// on `indent`: `(", ", ": ")` without, `(",", ": ")` with.
    pub separators: Option<(&'a str, &'a str)>,
    pub ensure_ascii: bool,
    pub sort_keys: bool,
    pub allow_nan: bool,
}

impl Default for Dump<'_> {
    /// `json.dumps(obj)` with no keywords.
    fn default() -> Self {
        Dump {
            indent: None,
            separators: None,
            ensure_ascii: true,
            sort_keys: false,
            allow_nan: true,
        }
    }
}

/// Why `dumps` or `loads` raised.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum JsonError {
    /// ValueError from `dumps(allow_nan=False)`.
    OutOfRangeFloat,
    /// JSONDecodeError: what was wrong and the byte offset.
    Decode(&'static str, usize),
    /// RecursionError: arrays and objects nested deeper than [`MAX_DEPTH`].
    TooDeep,
}

/// How deep `loads` nests arrays and objects before it gives up, about
/// where CPython's default recursion limit (1000) stops `json.loads`.
pub const MAX_DEPTH: usize = 1000;

impl JsonError {
    /// The Python exception class.
    pub fn exception(&self) -> &'static str {
        match self {
            JsonError::OutOfRangeFloat => "ValueError",
            JsonError::Decode(..) => "JSONDecodeError",
            JsonError::TooDeep => "RecursionError",
        }
    }
}

impl fmt::Display for JsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            JsonError::OutOfRangeFloat => {
                f.write_str("Out of range float values are not JSON compliant")
            }
            JsonError::Decode(what, at) => write!(f, "{what}: byte {at}"),
            JsonError::TooDeep => write!(f, "nested deeper than {MAX_DEPTH} levels"),
        }
    }
}

impl std::error::Error for JsonError {}

/// `json.dumps(value, **opts)`.
pub fn dumps(value: &Value, opts: &Dump<'_>) -> Result<String, JsonError> {
    let mut out = String::new();
    let (item, key) = opts.separators.unwrap_or(if opts.indent.is_some() {
        (",", ": ")
    } else {
        (", ", ": ")
    });
    let w = Writer { opts, item, key };
    w.value(&mut out, value, 0)?;
    Ok(out)
}

struct Writer<'o, 'a> {
    opts: &'o Dump<'a>,
    item: &'o str,
    key: &'o str,
}

impl Writer<'_, '_> {
    fn newline(&self, out: &mut String, level: usize) {
        if let Some(n) = self.opts.indent {
            out.push('\n');
            out.extend(std::iter::repeat_n(' ', n * level));
        }
    }

    /// The separator between two items, and the line break after it.
    fn between(&self, out: &mut String, level: usize) {
        out.push_str(self.item);
        self.newline(out, level);
    }

    fn value(&self, out: &mut String, v: &Value, level: usize) -> Result<(), JsonError> {
        match v {
            Value::Null => out.push_str("null"),
            Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Value::Int(i) => write!(out, "{i}").expect("write to String"),
            Value::Float(x) => out.push_str(&self.float(*x)?),
            Value::Str(s) => quote(out, s, self.opts.ensure_ascii),
            Value::Array(a) if a.is_empty() => out.push_str("[]"),
            Value::Array(a) => {
                out.push('[');
                self.newline(out, level + 1);
                for (i, x) in a.iter().enumerate() {
                    if i > 0 {
                        self.between(out, level + 1);
                    }
                    self.value(out, x, level + 1)?;
                }
                self.newline(out, level);
                out.push(']');
            }
            Value::Object(m) if m.is_empty() => out.push_str("{}"),
            Value::Object(m) => {
                let mut items: Vec<(&String, &Value)> = m.iter().collect();
                if self.opts.sort_keys {
                    items.sort_by(|a, b| a.0.cmp(b.0));
                }
                out.push('{');
                self.newline(out, level + 1);
                for (i, (k, x)) in items.into_iter().enumerate() {
                    if i > 0 {
                        self.between(out, level + 1);
                    }
                    quote(out, k, self.opts.ensure_ascii);
                    out.push_str(self.key);
                    self.value(out, x, level + 1)?;
                }
                self.newline(out, level);
                out.push('}');
            }
        }
        Ok(())
    }

    fn float(&self, x: f64) -> Result<String, JsonError> {
        if x.is_finite() {
            return Ok(repr(x));
        }
        if !self.opts.allow_nan {
            return Err(JsonError::OutOfRangeFloat);
        }
        Ok(if x.is_nan() {
            "NaN"
        } else if x > 0.0 {
            "Infinity"
        } else {
            "-Infinity"
        }
        .to_owned())
    }
}

/// A JSON string literal as Python's encoder writes it.
pub fn quote(out: &mut String, s: &str, ensure_ascii: bool) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\x08' => out.push_str("\\b"),
            '\x0c' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => write!(out, "\\u{:04x}", c as u32).expect("write to String"),
            c if ensure_ascii && !(' '..='~').contains(&c) => {
                let mut buf = [0u16; 2];
                for unit in c.encode_utf16(&mut buf) {
                    write!(out, "\\u{unit:04x}").expect("write to String");
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// `json.loads(text)`.
pub fn loads(text: &str) -> Result<Value, JsonError> {
    let mut p = Parser {
        s: text.as_bytes(),
        text,
        i: 0,
        depth: 0,
    };
    p.ws();
    let v = p.value()?;
    p.ws();
    if p.i != p.s.len() {
        return Err(JsonError::Decode("Extra data", p.i));
    }
    Ok(v)
}

struct Parser<'a> {
    s: &'a [u8],
    text: &'a str,
    i: usize,
    /// Arrays and objects open around `i`: a bound on the recursion.
    depth: usize,
}

impl Parser<'_> {
    fn err<T>(&self, what: &'static str) -> Result<T, JsonError> {
        Err(JsonError::Decode(what, self.i))
    }

    fn ws(&mut self) {
        while matches!(self.s.get(self.i), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.i += 1;
        }
    }

    fn eat(&mut self, lit: &str) -> bool {
        if self.s[self.i..].starts_with(lit.as_bytes()) {
            self.i += lit.len();
            true
        } else {
            false
        }
    }

    fn value(&mut self) -> Result<Value, JsonError> {
        match self.s.get(self.i) {
            Some(b'{' | b'[') => {
                if self.depth == MAX_DEPTH {
                    return Err(JsonError::TooDeep);
                }
                self.depth += 1;
                let v = if self.s[self.i] == b'{' {
                    self.object()
                } else {
                    self.array()
                };
                self.depth -= 1;
                v
            }
            Some(b'"') => Ok(Value::Str(self.string()?)),
            Some(b'n') if self.eat("null") => Ok(Value::Null),
            Some(b't') if self.eat("true") => Ok(Value::Bool(true)),
            Some(b'f') if self.eat("false") => Ok(Value::Bool(false)),
            Some(b'N') if self.eat("NaN") => Ok(Value::Float(f64::NAN)),
            Some(b'I') if self.eat("Infinity") => Ok(Value::Float(f64::INFINITY)),
            Some(b'-') if self.eat("-Infinity") => Ok(Value::Float(f64::NEG_INFINITY)),
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => self.err("Expecting value"),
        }
    }

    fn digits(&mut self) -> usize {
        let n = self.s[self.i..]
            .iter()
            .take_while(|c| c.is_ascii_digit())
            .count();
        self.i += n;
        n
    }

    /// `-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?`, as Python's scanner.
    fn number(&mut self) -> Result<Value, JsonError> {
        let start = self.i;
        if self.s[self.i] == b'-' {
            self.i += 1;
        }
        match self.s.get(self.i) {
            Some(b'0') => self.i += 1,
            Some(b'1'..=b'9') => {
                self.digits();
            }
            _ => {
                self.i = start;
                return self.err("Expecting value");
            }
        }
        let mut float = false;
        if self.s.get(self.i) == Some(&b'.')
            && self.s.get(self.i + 1).is_some_and(u8::is_ascii_digit)
        {
            self.i += 1;
            self.digits();
            float = true;
        }
        if matches!(self.s.get(self.i), Some(b'e' | b'E')) {
            let save = self.i;
            self.i += 1;
            if matches!(self.s.get(self.i), Some(b'+' | b'-')) {
                self.i += 1;
            }
            if self.digits() == 0 {
                self.i = save;
            } else {
                float = true;
            }
        }
        let lit = &self.text[start..self.i];
        if float {
            Ok(Value::Float(lit.parse().expect("a JSON number parses")))
        } else {
            lit.parse()
                .map(Value::Int)
                .or_else(|_| self.err("Integer beyond 128 bits"))
        }
    }

    fn hex4(&mut self) -> Result<u32, JsonError> {
        let h = self
            .s
            .get(self.i..self.i + 4)
            .and_then(|h| std::str::from_utf8(h).ok());
        match h.and_then(|h| {
            (h.bytes().all(|c| c.is_ascii_hexdigit())).then(|| u32::from_str_radix(h, 16).ok())
        }) {
            Some(Some(v)) => {
                self.i += 4;
                Ok(v)
            }
            _ => self.err("Invalid \\uXXXX escape"),
        }
    }

    fn string(&mut self) -> Result<String, JsonError> {
        self.i += 1;
        let mut out = String::new();
        loop {
            let run = self.s[self.i..]
                .iter()
                .take_while(|&&c| c != b'"' && c != b'\\' && c >= 0x20)
                .count();
            out.push_str(&self.text[self.i..self.i + run]);
            self.i += run;
            match self.s.get(self.i) {
                Some(b'"') => {
                    self.i += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    self.i += 1;
                    let c = match self.s.get(self.i) {
                        Some(b'"') => '"',
                        Some(b'\\') => '\\',
                        Some(b'/') => '/',
                        Some(b'b') => '\x08',
                        Some(b'f') => '\x0c',
                        Some(b'n') => '\n',
                        Some(b'r') => '\r',
                        Some(b't') => '\t',
                        Some(b'u') => {
                            self.i += 1;
                            let hi = self.hex4()?;
                            let cp = if (0xd800..0xdc00).contains(&hi)
                                && self.s[self.i..].starts_with(b"\\u")
                            {
                                let save = self.i;
                                self.i += 2;
                                let lo = self.hex4()?;
                                if (0xdc00..0xe000).contains(&lo) {
                                    0x10000 + ((hi - 0xd800) << 10) + (lo - 0xdc00)
                                } else {
                                    self.i = save;
                                    hi
                                }
                            } else {
                                hi
                            };
                            match char::from_u32(cp) {
                                Some(c) => {
                                    out.push(c);
                                    continue;
                                }
                                None => {
                                    return self
                                        .err("Lone surrogate (Python keeps it; this port cannot)");
                                }
                            }
                        }
                        _ => return self.err("Invalid \\escape"),
                    };
                    out.push(c);
                    self.i += 1;
                }
                Some(_) => return self.err("Invalid control character"),
                None => return self.err("Unterminated string"),
            }
        }
    }

    fn array(&mut self) -> Result<Value, JsonError> {
        self.i += 1;
        let mut a = Vec::new();
        self.ws();
        if self.eat("]") {
            return Ok(Value::Array(a));
        }
        loop {
            self.ws();
            a.push(self.value()?);
            self.ws();
            if self.eat("]") {
                return Ok(Value::Array(a));
            }
            if !self.eat(",") {
                return self.err("Expecting ',' delimiter");
            }
        }
    }

    fn object(&mut self) -> Result<Value, JsonError> {
        self.i += 1;
        let mut m = IndexMap::new();
        self.ws();
        if self.eat("}") {
            return Ok(Value::Object(m));
        }
        loop {
            self.ws();
            if self.s.get(self.i) != Some(&b'"') {
                return self.err("Expecting property name enclosed in double quotes");
            }
            let k = self.string()?;
            self.ws();
            if !self.eat(":") {
                return self.err("Expecting ':' delimiter");
            }
            self.ws();
            let v = self.value()?;
            // A duplicate key keeps its first place and takes the last value, as a dict does.
            m.insert(k, v);
            self.ws();
            if self.eat("}") {
                return Ok(Value::Object(m));
            }
            if !self.eat(",") {
                return self.err("Expecting ',' delimiter");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trip(text: &str, opts: &Dump<'_>) -> String {
        dumps(&loads(text).unwrap(), opts).unwrap()
    }

    #[test]
    fn dumps_matches_python() {
        let d = Dump::default();
        assert_eq!(
            round_trip(r#"{"a": [1, 2.0, "é😀\u007f"], "b": {}}"#, &d),
            r#"{"a": [1, 2.0, "\u00e9\ud83d\ude00\u007f"], "b": {}}"#
        );
        let one = Dump {
            indent: Some(1),
            ..Dump::default()
        };
        assert_eq!(
            round_trip(r#"{"b": [1, []], "a": null}"#, &one),
            "{\n \"b\": [\n  1,\n  []\n ],\n \"a\": null\n}"
        );
        let zero = Dump {
            indent: Some(0),
            ..Dump::default()
        };
        assert_eq!(round_trip("[1,2]", &zero), "[\n1,\n2\n]");
        let compact = Dump {
            separators: Some((",", ":")),
            sort_keys: true,
            ensure_ascii: false,
            ..Dump::default()
        };
        assert_eq!(
            round_trip(r#"{"b":1,"a":"é\u0001"}"#, &compact),
            "{\"a\":\"é\\u0001\",\"b\":1}"
        );
        assert_eq!(
            round_trip("[NaN, Infinity, -Infinity, 1e400, -0.0, 1E2]", &d),
            "[NaN, Infinity, -Infinity, Infinity, -0.0, 100.0]"
        );
        let strict = Dump {
            allow_nan: false,
            ..Dump::default()
        };
        assert_eq!(
            dumps(&Value::Float(f64::NAN), &strict),
            Err(JsonError::OutOfRangeFloat)
        );
    }

    #[test]
    fn loads_keeps_python_dict_semantics() {
        assert_eq!(
            round_trip(r#"{"a": 1, "b": 2, "a": 3}"#, &Dump::default()),
            r#"{"a": 3, "b": 2}"#
        );
        assert_eq!(loads("-0").unwrap(), Value::Int(0));
        assert_eq!(
            loads(" \"\\ud83d\\ude00\" ").unwrap(),
            Value::Str("😀".into())
        );
        for bad in [
            "",
            "[1,]",
            "{\"a\" 1}",
            "01",
            "1.",
            "\"a\u{1}\"",
            "[1] x",
            "{'a': 1}",
            "\"\\ud800\"",
        ] {
            assert!(loads(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn loads_stops_deep_nesting_with_an_error() {
        // 8 MiB, the main-thread stack og runs on: a test thread has 2 MiB,
        // too little for MAX_DEPTH levels in an unoptimised build.
        let t = std::thread::Builder::new().stack_size(8 << 20).spawn(|| {
            let nested = |n: usize| format!("{}{}", "[".repeat(n), "]".repeat(n));
            assert!(loads(&nested(MAX_DEPTH)).is_ok());
            assert_eq!(loads(&nested(MAX_DEPTH + 1)), Err(JsonError::TooDeep));
            // Far past the limit: an error, never a stack overflow.
            assert_eq!(loads(&"[".repeat(1_000_000)), Err(JsonError::TooDeep));
            assert_eq!(loads(&"{\"a\":".repeat(100_000)), Err(JsonError::TooDeep));
        });
        t.unwrap().join().unwrap();
        assert_eq!(JsonError::TooDeep.exception(), "RecursionError");
    }
}
