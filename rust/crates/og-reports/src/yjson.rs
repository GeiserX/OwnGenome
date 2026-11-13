//! A pull reader of JSON text: one event per call, no document tree.
//!
//! It takes what Python's `json.loads` takes (`NaN`, `Infinity` and
//! `-Infinity` included) and refuses what it refuses, with the byte offset.
//! The one difference: a lone UTF-16 surrogate in a string escape, which
//! Python keeps as a lone surrogate and Rust text cannot hold, is refused.

use std::borrow::Cow;

/// One step through the text.
#[derive(Debug, PartialEq)]
pub(crate) enum Event<'a> {
    StartObject,
    EndObject,
    StartArray,
    EndArray,
    /// An object key.
    Key(Cow<'a, str>),
    Str(Cow<'a, str>),
    /// A number, as written.
    Num(&'a str),
    Bool(bool),
    Null,
    /// `NaN`, `Infinity` or `-Infinity`.
    NonFinite,
    /// The end of the text, after the one top-level value.
    End,
}

impl Event<'_> {
    /// Python's `bool(v)` is false for this scalar: `null`, `false`, a zero
    /// number or an empty string.
    pub(crate) fn falsy(&self) -> bool {
        match self {
            Event::Null | Event::Bool(false) => true,
            Event::Num(n) => n.parse::<f64>().is_ok_and(|x| x == 0.0),
            Event::Str(s) => s.is_empty(),
            _ => false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Ctx {
    ObjFirst,
    ObjNext,
    ObjValue,
    ArrFirst,
    ArrNext,
}

/// The reader; [`Reader::next`] gives the events in text order.
pub(crate) struct Reader<'a> {
    s: &'a str,
    i: usize,
    stack: Vec<Ctx>,
    done: bool,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(s: &'a str) -> Self {
        Reader {
            s,
            i: 0,
            stack: Vec::new(),
            done: false,
        }
    }

    fn err<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("not valid JSON at byte {}: {what}", self.i))
    }

    fn ws(&mut self) {
        let b = self.s.as_bytes();
        while self.i < b.len() && matches!(b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.s.as_bytes().get(self.i).copied()
    }

    /// The next event, or why the text is not JSON.
    pub(crate) fn next(&mut self) -> Result<Event<'a>, String> {
        self.ws();
        let Some(&top) = self.stack.last() else {
            if self.done {
                return if self.i == self.s.len() {
                    Ok(Event::End)
                } else {
                    self.err("extra data")
                };
            }
            return self.value();
        };
        let n = self.stack.len() - 1;
        match top {
            Ctx::ObjFirst | Ctx::ObjNext => {
                if self.peek() == Some(b'}') {
                    self.i += 1;
                    return Ok(self.close(Event::EndObject));
                }
                if matches!(top, Ctx::ObjNext) {
                    if self.peek() != Some(b',') {
                        return self.err("expected ',' or '}'");
                    }
                    self.i += 1;
                    self.ws();
                }
                if self.peek() != Some(b'"') {
                    return self.err("expected a key in double quotes");
                }
                let key = self.string()?;
                self.ws();
                if self.peek() != Some(b':') {
                    return self.err("expected ':'");
                }
                self.i += 1;
                self.stack[n] = Ctx::ObjValue;
                Ok(Event::Key(key))
            }
            Ctx::ObjValue => {
                self.stack[n] = Ctx::ObjNext;
                self.value()
            }
            Ctx::ArrFirst | Ctx::ArrNext => {
                if self.peek() == Some(b']') {
                    self.i += 1;
                    return Ok(self.close(Event::EndArray));
                }
                if matches!(top, Ctx::ArrNext) {
                    if self.peek() != Some(b',') {
                        return self.err("expected ',' or ']'");
                    }
                    self.i += 1;
                    self.ws();
                }
                self.stack[n] = Ctx::ArrNext;
                self.value()
            }
        }
    }

    fn close(&mut self, e: Event<'a>) -> Event<'a> {
        self.stack.pop();
        if self.stack.is_empty() {
            self.done = true;
        }
        e
    }

    /// A scalar ends the text when it is the top-level value.
    fn scalar(&mut self, e: Event<'a>) -> Event<'a> {
        if self.stack.is_empty() {
            self.done = true;
        }
        e
    }

    fn literal(&mut self, word: &str, e: Event<'a>) -> Result<Event<'a>, String> {
        if self.s[self.i..].starts_with(word) {
            self.i += word.len();
            Ok(self.scalar(e))
        } else {
            self.err("expecting value")
        }
    }

    fn value(&mut self) -> Result<Event<'a>, String> {
        match self.peek() {
            Some(b'{') => {
                self.i += 1;
                self.stack.push(Ctx::ObjFirst);
                Ok(Event::StartObject)
            }
            Some(b'[') => {
                self.i += 1;
                self.stack.push(Ctx::ArrFirst);
                Ok(Event::StartArray)
            }
            Some(b'"') => {
                let s = self.string()?;
                Ok(self.scalar(Event::Str(s)))
            }
            Some(b't') => self.literal("true", Event::Bool(true)),
            Some(b'f') => self.literal("false", Event::Bool(false)),
            Some(b'n') => self.literal("null", Event::Null),
            Some(b'N') => self.literal("NaN", Event::NonFinite),
            Some(b'I') => self.literal("Infinity", Event::NonFinite),
            Some(b'-') if self.s[self.i..].starts_with("-Infinity") => {
                self.literal("-Infinity", Event::NonFinite)
            }
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => self.err("expecting value"),
        }
    }

    /// `-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][-+]?[0-9]+)?`, as Python's scanner
    /// reads it: a fraction or exponent without digits ends the number there.
    fn number(&mut self) -> Result<Event<'a>, String> {
        let b = self.s.as_bytes();
        let start = self.i;
        let digits = |i: &mut usize| {
            let s = *i;
            while *i < b.len() && b[*i].is_ascii_digit() {
                *i += 1;
            }
            *i - s
        };
        let mut i = self.i;
        if b[i] == b'-' {
            i += 1;
        }
        match b.get(i) {
            Some(b'0') => i += 1,
            Some(b'1'..=b'9') => {
                digits(&mut i);
            }
            _ => return self.err("expecting value"),
        }
        if b.get(i) == Some(&b'.') {
            let mut j = i + 1;
            if digits(&mut j) > 0 {
                i = j;
            }
        }
        if matches!(b.get(i), Some(b'e' | b'E')) {
            let mut j = i + 1;
            if matches!(b.get(j), Some(b'+' | b'-')) {
                j += 1;
            }
            if digits(&mut j) > 0 {
                i = j;
            }
        }
        self.i = i;
        Ok(self.scalar(Event::Num(&self.s[start..i])))
    }

    fn hex4(&mut self) -> Result<u32, String> {
        let h = self.s.get(self.i..self.i + 4).unwrap_or("");
        if h.len() != 4 || !h.bytes().all(|c| c.is_ascii_hexdigit()) {
            return self.err("invalid \\uXXXX escape");
        }
        self.i += 4;
        Ok(u32::from_str_radix(h, 16).unwrap_or(0))
    }

    /// A string at `"`; borrowed when it holds no escape.
    fn string(&mut self) -> Result<Cow<'a, str>, String> {
        let b = self.s.as_bytes();
        self.i += 1;
        let start = self.i;
        let mut owned: Option<String> = None;
        let mut run = start;
        loop {
            let Some(&c) = b.get(self.i) else {
                return self.err("unterminated string");
            };
            match c {
                b'"' => {
                    let tail = &self.s[run..self.i];
                    self.i += 1;
                    return Ok(match owned {
                        None => Cow::Borrowed(tail),
                        Some(mut o) => {
                            o.push_str(tail);
                            Cow::Owned(o)
                        }
                    });
                }
                b'\\' => {
                    let o = owned.get_or_insert_with(String::new);
                    o.push_str(&self.s[run..self.i]);
                    self.i += 1;
                    let Some(&e) = b.get(self.i) else {
                        return self.err("unterminated string");
                    };
                    self.i += 1;
                    let ch = match e {
                        b'"' => '"',
                        b'\\' => '\\',
                        b'/' => '/',
                        b'b' => '\u{8}',
                        b'f' => '\u{c}',
                        b'n' => '\n',
                        b'r' => '\r',
                        b't' => '\t',
                        b'u' => {
                            let hi = self.hex4()?;
                            let cp = if (0xD800..0xDC00).contains(&hi)
                                && self.s[self.i..].starts_with("\\u")
                            {
                                let save = self.i;
                                self.i += 2;
                                let lo = self.hex4()?;
                                if (0xDC00..0xE000).contains(&lo) {
                                    0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)
                                } else {
                                    self.i = save;
                                    hi
                                }
                            } else {
                                hi
                            };
                            match char::from_u32(cp) {
                                Some(ch) => ch,
                                None => return self.err("a lone surrogate in a \\u escape"),
                            }
                        }
                        _ => return self.err("invalid escape"),
                    };
                    let o = owned.get_or_insert_with(String::new);
                    o.push(ch);
                    run = self.i;
                }
                0..=0x1f => return self.err("invalid control character in a string"),
                _ => self.i += 1,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(s: &str) -> Result<Vec<String>, String> {
        let mut r = Reader::new(s);
        let mut out = Vec::new();
        loop {
            let e = r.next()?;
            if e == Event::End {
                return Ok(out);
            }
            out.push(format!("{e:?}"));
        }
    }

    #[test]
    fn events_in_order() {
        let got = all(
            r#" {"a": [1, -2.5e3, true, null, NaN, -Infinity], "b\u00e9\n": "x\"y", "c": {}} "#,
        )
        .unwrap();
        assert_eq!(
            got,
            [
                "StartObject",
                "Key(\"a\")",
                "StartArray",
                "Num(\"1\")",
                "Num(\"-2.5e3\")",
                "Bool(true)",
                "Null",
                "NonFinite",
                "NonFinite",
                "EndArray",
                "Key(\"bé\\n\")",
                "Str(\"x\\\"y\")",
                "Key(\"c\")",
                "StartObject",
                "EndObject",
                "EndObject"
            ]
        );
        assert_eq!(all("\"\\ud83d\\ude00\"").unwrap(), ["Str(\"😀\")"]);
        assert_eq!(all("0").unwrap(), ["Num(\"0\")"]);
    }

    #[test]
    fn refusals() {
        for bad in [
            "",
            "{",
            "[1,]",
            "{\"a\" 1}",
            "{\"a\":1,}",
            "[1 2]",
            "01",
            "{} x",
            "\"a\nb\"",
            "'a'",
            "\"\\ud800\"",
            "\"\\x\"",
            "tru",
            "[-]",
            "{1:2}",
        ] {
            assert!(all(bad).is_err(), "{bad:?} should be refused");
        }
    }

    #[test]
    fn falsy_as_python() {
        for (e, f) in [
            (Event::Null, true),
            (Event::Bool(false), true),
            (Event::Num("0"), true),
            (Event::Num("-0.0e5"), true),
            (Event::Num("1"), false),
            (Event::Str(Cow::Borrowed("")), true),
            (Event::Str(Cow::Borrowed("x")), false),
            (Event::NonFinite, false),
        ] {
            assert_eq!(e.falsy(), f, "{e:?}");
        }
    }
}
