//! Image rules keep their original linear matcher. Personal-data rules use
//! a Python-subset compiler and a backtracking matcher with assertions.
//!
//! Every look-around of the image scanner is a single character class: "not
//! after one of these" before the match and "not followed by one of these"
//! after it. A [`Rule`] keeps the matched text in capture group 1 and
//! writes the look-ahead as a consumed `(?:[^...]|$)` after it, so the
//! engine still finds, at each start, the match Python's backtracking would
//! pick. The look-behind is the character before group 1, checked in code;
//! when it fails, the search starts again one character later, as Python
//! tries the next start. The next match is searched from the end of group
//! 1, so the character the look-ahead read is not used up.

use regex::Regex;

/// A pattern with its look-behind.
pub struct Rule {
    re: Regex,
    /// True for a character the match must not follow.
    not_after: Option<fn(char) -> bool>,
}

impl Rule {
    /// `pattern` holds the match in group 1, or is the match when it has no
    /// group.
    pub fn new(pattern: &str, not_after: Option<fn(char) -> bool>) -> Result<Rule, regex::Error> {
        Ok(Rule {
            re: Regex::new(pattern)?,
            not_after,
        })
    }

    /// The spans of the image rule over `line`, in order.
    pub fn find_iter<'a>(&'a self, line: &'a str) -> impl Iterator<Item = (usize, usize)> + 'a {
        let mut pos = 0;
        std::iter::from_fn(move || {
            while pos <= line.len() {
                let caps = self.re.captures_at(line, pos)?;
                let m = caps.get(1).or_else(|| caps.get(0))?;
                let (start, end) = (m.start(), m.end());
                let refused = self
                    .not_after
                    .is_some_and(|bad| line[..start].chars().next_back().is_some_and(bad));
                if refused {
                    pos = next_char(line, start);
                    continue;
                }
                pos = if end > start {
                    end
                } else {
                    next_char(line, end)
                };
                return Some((start, end));
            }
            None
        })
    }
}

/// The byte offset of the character after `i`, or one past the end.
fn next_char(s: &str, i: usize) -> usize {
    s[i..]
        .chars()
        .next()
        .map_or(s.len() + 1, |c| i + c.len_utf8())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn digit_or_dot(c: char) -> bool {
        c.is_ascii_digit() || c == '.'
    }

    #[test]
    fn behind_and_ahead_as_python() {
        // Python: re.finditer(r'(?<![0-9.])[0-9]{1,4}/[0-9]{1,4}(?![0-9/])', s)
        let r = Rule::new(r"([0-9]{1,4}/[0-9]{1,4})(?:[^0-9/]|$)", Some(digit_or_dot)).unwrap();
        fn spans<'a>(r: &Rule, s: &'a str) -> Vec<&'a str> {
            r.find_iter(s).map(|(a, b)| &s[a..b]).collect()
        }
        assert_eq!(spans(&r, "17/19"), ["17/19"]);
        assert_eq!(spans(&r, "1.17/19 x"), Vec::<&str>::new());
        assert_eq!(spans(&r, "17/19/3"), ["19/3"]);
        assert_eq!(spans(&r, "a 17/19 and 2/3"), ["17/19", "2/3"]);
        // The look-ahead character is not used up: two matches side by side.
        let w = Rule::new(r"(ab)(?:[^b]|$)", None).unwrap();
        assert_eq!(w.find_iter("abab").count(), 2);
    }

    #[test]
    fn a_shorter_match_where_the_longest_fails_ahead() {
        // Python finds "a/b:1" in "a/b:1.2/x": the tag backs off to "1",
        // which a "." follows.
        let r = Rule::new(r"([a-z]+/[a-z]+:[0-9][0-9.]*)(?:[^0-9a-z/]|$)", None).unwrap();
        let s = "a/b:1.2/x";
        let got: Vec<_> = r.find_iter(s).map(|(a, b)| &s[a..b]).collect();
        assert_eq!(got, ["a/b:1"]);
    }
}

// Personal-data patterns have one translation and one search path. Assertions
// stay in the engine so a failed assertion can backtrack at the same start.
const WORD: &str = r"[\p{L}\p{N}_]";
const SPACE: &str = r"[\x09-\x0d\x1c-\x20\x85\x{a0}\x{1680}\x{2000}-\x{200a}\x{2028}\x{2029}\x{202f}\x{205f}\x{3000}]";

/// Why a pattern was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Refusal {
    /// It names, under an active `i` flag, a letter whose case pairing
    /// differs between Unicode versions.
    CaseDependent,
    Other,
}

#[derive(Clone)]
pub struct PythonRule {
    re: fancy_regex::Regex,
    after_empty: fancy_regex::Regex,
}

impl PythonRule {
    pub fn compile(pattern: &str, builtin: bool) -> Result<Self, ()> {
        Self::compile_explained(pattern, builtin).map_err(|_| ())
    }

    /// As `compile`, and says when the refusal is the case-pairing guard, so
    /// the caller can tell the user what to change.
    pub fn compile_explained(pattern: &str, builtin: bool) -> Result<Self, Refusal> {
        Self::compile_inner(pattern, builtin).map_err(|case| {
            if case {
                Refusal::CaseDependent
            } else {
                Refusal::Other
            }
        })
    }

    /// The error is true when the case-pairing guard refused the pattern.
    fn compile_inner(pattern: &str, builtin: bool) -> Result<Self, bool> {
        if pattern.chars().any(crate::personal_data::version_dependent) {
            return Err(false);
        }
        if pattern.len() > 16_384 {
            return Err(false);
        }
        let mut parser = PythonPattern {
            chars: pattern.chars().collect(),
            at: 0,
            builtin,
            depth: 0,
            names: std::collections::HashSet::new(),
            case_refused: false,
        };
        let (translated, _) = parser
            .group(Flags::default(), false)
            .map_err(|()| parser.case_refused)?;
        // \G is the search offset. After yielding an empty match, forbid only
        // a match ending at that offset; the engine then tries nonempty
        // alternatives at the same start, exactly as Python finditer does.
        let after_empty = format!(r"(?:{translated})(?!\G)");
        let build = |p: &str| {
            fancy_regex::RegexBuilder::new(p)
                .backtrack_limit(1_000_000)
                .build()
                .map_err(|_| false)
        };
        Ok(Self {
            re: build(&translated)?,
            after_empty: build(&after_empty)?,
        })
    }

    pub fn find_iter<'a>(
        &'a self,
        text: &'a str,
    ) -> impl Iterator<Item = Result<(usize, usize), ()>> + 'a {
        let (mut pos, mut empty, mut done) = (0, false, false);
        std::iter::from_fn(move || {
            if done {
                return None;
            }
            let re = if empty { &self.after_empty } else { &self.re };
            match re.find_from_pos(text, pos) {
                Ok(Some(m)) => {
                    pos = m.end();
                    empty = m.start() == m.end();
                    Some(Ok((m.start(), m.end())))
                }
                Ok(None) => {
                    done = true;
                    None
                }
                Err(_) => {
                    done = true;
                    Some(Err(()))
                }
            }
        })
    }

    pub fn is_match(&self, text: &str) -> Result<bool, ()> {
        self.find_iter(text).next().transpose().map(|m| m.is_some())
    }
}

#[derive(Clone, Copy, Default)]
struct Flags {
    i: bool,
    m: bool,
    s: bool,
    x: bool,
}

struct PythonPattern {
    chars: Vec<char>,
    at: usize,
    builtin: bool,
    depth: usize,
    names: std::collections::HashSet<String>,
    /// Set when a refusal comes from the case-pairing guard.
    case_refused: bool,
}

impl PythonPattern {
    /// `literal`, recording a refusal by the case-pairing guard.
    fn lit(&mut self, c: char, ci: bool) -> Result<String, ()> {
        self.case_refused |= ci && crate::personal_data::range_regex_case_dependent(c, c);
        literal(c, ci)
    }

    fn take(&mut self) -> Result<char, ()> {
        let c = self.chars.get(self.at).copied().ok_or(())?;
        self.at += 1;
        Ok(c)
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.at).copied()
    }

    fn group(&mut self, mut flags: Flags, nested: bool) -> Result<(String, bool), ()> {
        let mut out = String::new();
        // Nullability is structural: concatenate with AND, alternate with OR.
        // Keep the last atom separate until its quantifier is known.
        let (mut prefix_nullable, mut atom_nullable, mut alt_nullable) = (true, true, false);
        let (mut repeatable, mut lazy) = (false, false);
        while let Some(c) = self.peek() {
            self.at += 1;
            if flags.x && matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{b}' | '\u{c}') {
                continue;
            }
            if flags.x && c == '#' {
                while self.peek().is_some_and(|c| c != '\n') {
                    self.at += 1;
                }
                continue;
            }
            if !matches!(c, ')' | '|' | '*' | '+' | '?' | '{') {
                prefix_nullable &= atom_nullable;
                atom_nullable = true;
            }
            match c {
                ')' if nested => {
                    return Ok((out, alt_nullable || (prefix_nullable && atom_nullable)));
                }
                ')' => return Err(()),
                '(' => {
                    repeatable = true;
                    lazy = false;
                    if self.depth >= 64 {
                        return Err(());
                    }
                    let mut prefix = "?:".to_owned();
                    let mut inner = flags;
                    if self.peek() == Some('?') {
                        self.at += 1;
                        match self.take()? {
                            ':' => (),
                            '=' if self.builtin => prefix = "?=".into(),
                            '!' if self.builtin => prefix = "?!".into(),
                            '<' if self.builtin => {
                                let c = self.take()?;
                                if !matches!(c, '=' | '!') {
                                    return Err(());
                                }
                                prefix = format!("?<{c}");
                            }
                            'P' => {
                                if self.take()? != '<' {
                                    return Err(());
                                }
                                let mut name = String::new();
                                loop {
                                    let c = self.take()?;
                                    if c == '>' {
                                        break;
                                    }
                                    name.push(c);
                                }
                                if name.is_empty()
                                    || !name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
                                    || name.starts_with(|c: char| c.is_ascii_digit())
                                    || !self.names.insert(name.clone())
                                {
                                    return Err(());
                                }
                            }
                            c if "imsx-".contains(c) => {
                                self.at -= 1;
                                let mut on = true;
                                let mut seen = std::collections::HashSet::new();
                                let mut disabled = false;
                                loop {
                                    let flag = self.take()?;
                                    if matches!(flag, 'i' | 'm' | 's' | 'x') {
                                        if !seen.insert(flag) {
                                            return Err(());
                                        }
                                        disabled |= !on;
                                    }
                                    match flag {
                                        '-' if on => on = false,
                                        'i' => inner.i = on,
                                        'm' => inner.m = on,
                                        's' => inner.s = on,
                                        'x' => inner.x = on,
                                        ':' if !seen.is_empty() && (on || disabled) => break,
                                        ')' if on
                                            && !seen.is_empty()
                                            && out.is_empty()
                                            && !nested =>
                                        {
                                            flags = inner;
                                            break;
                                        }
                                        _ => return Err(()),
                                    }
                                    if self.chars[self.at - 1] == ')' {
                                        break;
                                    }
                                }
                                if self.chars[self.at - 1] == ')' {
                                    repeatable = false;
                                    continue;
                                }
                            }
                            _ => return Err(()),
                        }
                    }
                    out.push('(');
                    out.push_str(&prefix);
                    self.depth += 1;
                    let (body, nullable) = self.group(inner, true)?;
                    atom_nullable = prefix != "?:" || nullable;
                    out.push_str(&body);
                    self.depth -= 1;
                    out.push(')');
                }
                '[' => {
                    out.push_str(&self.class(flags.i)?);
                    atom_nullable = false;
                    repeatable = true;
                    lazy = false;
                }
                '\\' => {
                    repeatable = !matches!(self.peek(), Some('b' | 'B' | 'A' | 'Z'));
                    atom_nullable = !repeatable;
                    lazy = false;
                    out.push_str(&self.escape(flags.i, false)?);
                }
                '.' => {
                    out.push_str(if flags.s { r"[\s\S]" } else { r"[^\n]" });
                    atom_nullable = false;
                    repeatable = true;
                    lazy = false;
                }
                '^' => {
                    out.push_str(if flags.m { r"(?:\A|(?<=\n))" } else { r"\A" });
                    repeatable = false;
                    lazy = false;
                }
                '$' => {
                    out.push_str(if flags.m { r"(?=\n|\z)" } else { r"(?=\n?\z)" });
                    repeatable = false;
                    lazy = false;
                }
                '|' => {
                    alt_nullable |= prefix_nullable && atom_nullable;
                    prefix_nullable = true;
                    atom_nullable = true;
                    out.push('|');
                    repeatable = false;
                    lazy = false;
                }
                '*' | '+' | '?' => {
                    if c == '?' && lazy {
                        out.push(c);
                        lazy = false;
                        continue;
                    }
                    if !repeatable || atom_nullable {
                        return Err(());
                    }
                    repeatable = false;
                    lazy = true;
                    // Possessive quantifiers are outside the supported subset.
                    if matches!(c, '*' | '+' | '?') && self.peek() == Some('+') {
                        return Err(());
                    }
                    atom_nullable = c != '+';
                    out.push(c);
                }
                '{' => {
                    if !repeatable || atom_nullable {
                        return Err(());
                    }
                    repeatable = false;
                    lazy = true;
                    out.push(c);
                    let minimum_at = self.at;
                    loop {
                        let c = self.take()?;
                        out.push(c);
                        if c == '}' {
                            break;
                        }
                        if !c.is_ascii_digit() && c != ',' {
                            return Err(());
                        }
                    }
                    // Python permits {,n}; a missing minimum means zero.
                    atom_nullable = self.chars[minimum_at..self.at - 1]
                        .iter()
                        .take_while(|&&c| c != ',')
                        .all(|&c| c == '0');
                    if self.peek() == Some('+') {
                        return Err(());
                    }
                }
                '}' => return Err(()),
                c => {
                    out.push_str(&self.lit(c, flags.i)?);
                    atom_nullable = false;
                    repeatable = true;
                    lazy = false;
                }
            }
        }
        if nested {
            Err(())
        } else {
            Ok((out, alt_nullable || (prefix_nullable && atom_nullable)))
        }
    }

    fn escape(&mut self, ci: bool, class: bool) -> Result<String, ()> {
        let c = self.take()?;
        Ok(match c {
            'w' => WORD.into(),
            'W' => r"[^\p{L}\p{N}_]".into(),
            's' => SPACE.into(),
            'S' => format!("[^{SPACE}]"),
            'd' => r"[\p{Nd}]".into(),
            'D' => r"[^\p{Nd}]".into(),
            'b' if class => r"\x08".into(),
            'b' => format!(r"(?:(?<!{WORD})(?={WORD})|(?<={WORD})(?!{WORD}))"),
            'B' if !class => format!(r"(?:(?<={WORD})(?={WORD})|(?<!{WORD})(?!{WORD}))"),
            'A' if !class => r"\A".into(),
            'Z' if !class => r"\z".into(),
            'a' => r"\x07".into(),
            'f' => r"\x0c".into(),
            'n' => r"\x0a".into(),
            'r' => r"\x0d".into(),
            't' => r"\x09".into(),
            'v' => r"\x0b".into(),
            'x' | 'u' | 'U' => {
                let n = match c {
                    'x' => 2,
                    'u' => 4,
                    _ => 8,
                };
                let mut hex = String::new();
                for _ in 0..n {
                    let c = self.take()?;
                    if !c.is_ascii_hexdigit() {
                        return Err(());
                    }
                    hex.push(c);
                }
                let c = char::from_u32(u32::from_str_radix(&hex, 16).map_err(|_| ())?).ok_or(())?;
                self.lit(c, ci)?
            }
            c if c.is_ascii_alphanumeric() => return Err(()),
            c => self.lit(c, ci)?,
        })
    }

    fn class(&mut self, ci: bool) -> Result<String, ()> {
        let negative = self.peek() == Some('^');
        if negative {
            self.at += 1;
        }
        let mut atoms = Vec::new();
        let mut first = true;
        loop {
            let c = self.take()?;
            if c == ']' && !first {
                break;
            }
            first = false;
            match c {
                '\\' => atoms.push((self.escape(false, true)?, false)),
                '[' => return Err(()),
                '-' => atoms.push((literal('-', false)?, true)),
                c => atoms.push((literal(c, false)?, false)),
            }
        }
        let mut body = String::new();
        let mut i = 0;
        while i < atoms.len() {
            let atom = &atoms[i].0;
            if i + 2 < atoms.len() && atoms[i + 1].1 {
                // Categories are not case folded in Python, even inside
                // classes. Only literal range endpoints can form a range.
                let start = class_literal(atom).ok_or(())?;
                let end = class_literal(&atoms[i + 2].0).ok_or(())?;
                self.case_refused |=
                    ci && crate::personal_data::range_regex_case_dependent(start, end);
                if start > end
                    || crate::personal_data::range_version_dependent(start, end)
                    || self.case_refused
                {
                    return Err(());
                }
                let range = format!("[{atom}-{}]", atoms[i + 2].0);
                body.push_str(&if ci { fold_class(&range)? } else { range });
                i += 3;
            } else {
                if atoms[i].1 && i != 0 && i + 1 != atoms.len() {
                    return Err(());
                }
                if ci && let Some(c) = class_literal(atom) {
                    if crate::personal_data::range_regex_case_dependent(c, c) {
                        self.case_refused = true;
                        return Err(());
                    }
                    body.push_str(&fold_class(&format!("[{atom}]"))?);
                } else {
                    body.push_str(atom);
                }
                i += 1;
            }
        }
        Ok(format!("[{}{body}]", if negative { "^" } else { "" }))
    }
}

fn class_literal(atom: &str) -> Option<char> {
    let hir = regex_syntax::Parser::new().parse(atom).ok()?;
    let regex_syntax::hir::HirKind::Literal(value) = hir.kind() else {
        return None;
    };
    let text = std::str::from_utf8(&value.0).ok()?;
    (text.chars().count() == 1)
        .then(|| text.chars().next())
        .flatten()
}

fn literal(c: char, ci: bool) -> Result<String, ()> {
    if crate::personal_data::version_dependent(c)
        || (ci && crate::personal_data::range_regex_case_dependent(c, c))
    {
        return Err(());
    }
    let escaped = format!(r"\x{{{:x}}}", u32::from(c));
    if ci {
        fold_class(&format!("[{escaped}]"))
    } else {
        Ok(escaped)
    }
}

fn fold_class(class: &str) -> Result<String, ()> {
    use regex_syntax::hir::{Class, HirKind};
    // Rust's simple fold follows Python except the extra dotted/dotless I.
    let hir = regex_syntax::ParserBuilder::new()
        .case_insensitive(true)
        .build()
        .parse(class)
        .map_err(|_| ())?;
    let set = match hir.kind() {
        HirKind::Class(Class::Unicode(set)) => set.clone(),
        HirKind::Literal(value) => {
            let text = std::str::from_utf8(&value.0).map_err(|_| ())?;
            let c = text.chars().next().ok_or(())?;
            if text.chars().count() != 1 {
                return Err(());
            }
            regex_syntax::hir::ClassUnicode::new([regex_syntax::hir::ClassUnicodeRange::new(c, c)])
        }
        _ => return Err(()),
    };
    let mut out = String::from("[");
    let mut i = false;
    for range in set.iter() {
        out.push_str(&format!(
            r"\x{{{:x}}}-\x{{{:x}}}",
            u32::from(range.start()),
            u32::from(range.end())
        ));
        i |= range.start() <= 'I' && 'I' <= range.end()
            || range.start() <= 'i' && 'i' <= range.end()
            || range.start() <= 'İ' && 'İ' <= range.end()
            || range.start() <= 'ı' && 'ı' <= range.end();
    }
    if i {
        out.push_str(r"Ii\x{130}\x{131}");
    }
    out.push(']');
    Ok(out)
}

#[cfg(test)]
mod python_tests {
    use super::*;

    #[test]
    fn finditer_spans_match_real_python() {
        for pattern in [
            r"^|a",
            r"a*?",
            r"a|",
            r"\w*",
            r"\b|\W+",
            r"a$",
            r"(?m)^|$",
            r"(?s:.)*?",
            r"(?i)[^I]+",
            r"(?i:İ|ı|S|k)+",
            r"[\s\S]",
            r"(?i)[^\W]+",
            r"(?i)[\w]+",
            r"\D*\d?",
            r"\Bz",
            r"(?x) a + # comment",
        ] {
            for text in ["", "aa", "a\n", "a\n\n", "\u{301}aİıſK9\u{1c}\u{345}Ιι"] {
                let output = std::process::Command::new("python3").args(["-I", "-c",
                    "import re,sys; print(','.join(str(m.start())+':'+str(m.end()) for m in re.finditer(sys.argv[1],sys.argv[2])))", pattern, text]).output().unwrap();
                assert!(output.status.success(), "span oracle failed");
                let rule = PythonRule::compile(pattern, false).expect("supported span pattern");
                let actual = rule
                    .find_iter(text)
                    .map(|m| {
                        let (a, b) = m.unwrap();
                        format!(
                            "{}:{}",
                            text[..a].chars().count(),
                            text[..b].chars().count()
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                assert_eq!(actual, String::from_utf8(output.stdout).unwrap().trim());
            }
        }
    }

    #[test]
    fn round3_nullable_atoms_are_refused() {
        for (pattern, text, spans) in [
            ("(a*|b*?){1,3}", "ba", "0:0,0:2,2:2"),
            ("(?:a*|[ab]{0,2}){0,3}", "baa", "0:0,0:3,3:3"),
        ] {
            let output = std::process::Command::new("python3").args(["-I", "-c",
                "import re,sys; print(','.join(str(m.start())+':'+str(m.end()) for m in re.finditer(sys.argv[1],sys.argv[2])))", pattern, text]).output().unwrap();
            assert!(output.status.success());
            assert_eq!(String::from_utf8(output.stdout).unwrap().trim(), spans);
            assert!(
                PythonRule::compile(pattern, false).is_err(),
                "unsupported span case accepted"
            );
        }
        for pattern in [
            "(a*|b*?){1,3}",
            "(?:a*|[ab]{0,2}){0,3}",
            "(?:a?)+",
            "(?:a*)?",
            "(?:|a)*?",
            "()*",
            "(?:a{0}){2}",
            "(?i:(?:a?|b)){1,2}",
            "(?:^|a)+",
            "(?:a|$)?",
        ] {
            assert!(
                PythonRule::compile(pattern, false).is_err(),
                "nullable operand accepted"
            );
        }
    }

    // A finite grammar, exhaustively enumerated by AST node count. Every
    // terminal costs one node; each unary/binary operator costs one more.
    // No random selection or accepted-pattern sampling occurs in this model.
    const MODEL_ATOMS: &[&str] = &[
        "", "a", "b", "0", " ", "\u{301}", ".", "[ab]", "[^ab]", "[a-b]", r"\w", r"\W", r"\d",
        r"\D", r"\s", r"\S", r"\b", r"\B", "^", "$", r"\A", r"\Z", r"\x61", r"\u0301", r"[\s\S]",
        r"[\w0]",
    ];
    const MODEL_GROUPS: &[&str] = &["", "?:", "?P<n>", "?i:", "?m:", "?s:", "?x:", "?-i:"];
    const MODEL_QUANTIFIERS: &[&str] = &[
        "?", "*", "+", "??", "*?", "+?", "{0}", "{1}", "{0,2}", "{1,3}", "{2,}", "{,2}", "{0,2}?",
        "{1,3}?",
    ];

    fn unary_patterns(p: &str) -> impl Iterator<Item = String> + '_ {
        MODEL_GROUPS
            .iter()
            .map(move |g| format!("({g}{p})"))
            .chain(MODEL_QUANTIFIERS.iter().map(move |q| format!("(?:{p}){q}")))
    }

    fn model_patterns(size: usize) -> Vec<String> {
        use std::collections::BTreeSet;
        let mut levels = vec![
            BTreeSet::new(),
            MODEL_ATOMS
                .iter()
                .map(|s| s.to_string())
                .collect::<BTreeSet<_>>(),
        ];
        for n in 2..=size {
            let mut level = BTreeSet::new();
            for p in &levels[n - 1] {
                level.extend(unary_patterns(p));
            }
            for left in 1..n - 1 {
                for a in &levels[left] {
                    for b in &levels[n - 1 - left] {
                        level.insert(format!("(?:{a})(?:{b})"));
                        level.insert(format!("(?:{a}|{b})"));
                    }
                }
            }
            levels.push(level);
        }
        levels
            .into_iter()
            .flatten()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    fn model_texts(length: usize) -> Vec<String> {
        let mut all = vec![String::new()];
        let mut level = vec![String::new()];
        for _ in 0..length {
            level = level
                .iter()
                .flat_map(|p| ['a', 'b', '0', ' ', '\u{301}', 'é'].map(|c| format!("{p}{c}")))
                .collect();
            all.extend(level.iter().cloned());
        }
        all
    }

    fn span_list(rule: &PythonRule, text: &str) -> String {
        rule.find_iter(text)
            .map(|m| {
                let (a, b) = m.expect("model matcher exhausted");
                format!(
                    "{}:{}",
                    text[..a].chars().count(),
                    text[..b].chars().count()
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    }

    fn compare_model(label: &str, patterns: &[String], texts: &[String]) -> usize {
        use std::io::Write;
        use std::process::{Command, Stdio};
        let (mut accepted, mut refused) = (0, 0);
        for (batch, chunk) in patterns.chunks(128).enumerate() {
            let rules: Vec<_> = chunk
                .iter()
                .filter_map(|p| match PythonRule::compile(p, false) {
                    Ok(r) => Some((p, r)),
                    Err(()) => {
                        refused += 1;
                        None
                    }
                })
                .collect();
            if rules.is_empty() {
                continue;
            }
            accepted += rules.len();
            let array = |items: Vec<&str>| {
                let mut out = String::from("[");
                for (i, item) in items.iter().enumerate() {
                    if i != 0 {
                        out.push(',');
                    }
                    pycompat::json::quote(&mut out, item, true);
                }
                out.push(']');
                out
            };
            let input = format!(
                "[{},{}]",
                array(rules.iter().map(|(p, _)| p.as_str()).collect()),
                array(texts.iter().map(String::as_str).collect())
            );
            let mut child = Command::new("python3")
                .args([
                    "-I",
                    "-c",
                    r#"
import json,re,sys
patterns,texts=json.load(sys.stdin)
for p in patterns:
    r=re.compile(p)
    for text in texts:
        print(','.join(str(m.start())+':'+str(m.end()) for m in r.finditer(text)))
"#,
                ])
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            let mut stdin = child.stdin.take().unwrap();
            let writer = std::thread::spawn(move || stdin.write_all(input.as_bytes()));
            let output = child.wait_with_output().unwrap();
            writer.join().unwrap().unwrap();
            assert!(output.status.success(), "model oracle failed");
            assert!(output.stderr.is_empty(), "model oracle diagnostic");
            let oracle = String::from_utf8(output.stdout).unwrap();
            let mut spans = oracle.lines();
            for (index, (_, rule)) in rules.iter().enumerate() {
                for (t, text) in texts.iter().enumerate() {
                    assert_eq!(
                        Some(span_list(rule, text).as_str()),
                        spans.next(),
                        "{label}: batch {batch}, accepted pattern {index}, text {t}"
                    );
                }
            }
            assert!(spans.next().is_none(), "extra oracle rows");
        }
        let pairs = accepted * texts.len();
        println!(
            "{label}: {} patterns, {accepted} accepted, {refused} refused, {} texts, {pairs} pairs, zero span differences",
            patterns.len(),
            texts.len()
        );
        pairs
    }

    #[test]
    fn exhaustive_small_model_spans() {
        // CI: every <=2-node pattern against all <=2-character texts.
        assert!(compare_model("exhaustive-ci", &model_patterns(2), &model_texts(2)) > 10_000);
    }

    #[test]
    #[ignore = "large local model; run explicitly with --ignored --nocapture"]
    fn exhaustive_large_model_spans() {
        // Local: every <=3-node pattern against all <=3-character texts.
        assert!(compare_model("exhaustive-local", &model_patterns(3), &model_texts(3)) > 1_000_000);
    }

    fn fuzz_pattern(state: &mut u64, depth: usize) -> String {
        fn next(state: &mut u64) -> usize {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            *state as usize
        }
        if depth == 0 || next(state).is_multiple_of(4) {
            return MODEL_ATOMS[next(state) % MODEL_ATOMS.len()].into();
        }
        let p = fuzz_pattern(state, depth - 1);
        match next(state) % 4 {
            0 => format!("({}{p})", MODEL_GROUPS[next(state) % MODEL_GROUPS.len()]),
            1 => format!(
                "(?:{p}){}",
                MODEL_QUANTIFIERS[next(state) % MODEL_QUANTIFIERS.len()]
            ),
            2 => format!("(?:{p}|{})", fuzz_pattern(state, depth - 1)),
            _ => format!("(?:{p})(?:{})", fuzz_pattern(state, depth - 1)),
        }
    }

    #[test]
    fn round3_two_seed_span_fuzz() {
        for mut seed in [0x401ad33, 0x3ad401] {
            let label = format!("fuzz-{seed:x}");
            let mut patterns = Vec::new();
            // Focus the same seeds on nullable alternatives and nested repeats.
            for _ in 0..1500 {
                let p = fuzz_pattern(&mut seed, 4);
                patterns.push(p.clone());
                patterns.push(format!("(?:{p}|a*?){{1,3}}"));
            }
            assert!(compare_model(&label, &patterns, &model_texts(3)) > 50_000);
        }
    }
}
