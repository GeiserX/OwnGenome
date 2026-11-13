//! Python's one-character look-behind and look-ahead on the `regex` crate,
//! which has neither.
//!
//! Every look-around of the scanners is a single character class: "not
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
    /// Capture describing the match span: 1 for rewritten look-around,
    /// 0 for a user's unmodified regular expression.
    capture: usize,
    /// True for a character the match must not follow.
    not_after: Option<fn(char) -> bool>,
}

impl Rule {
    /// `pattern` holds the match in group 1, or is the match when it has no
    /// group.
    pub fn new(pattern: &str, not_after: Option<fn(char) -> bool>) -> Result<Rule, regex::Error> {
        Ok(Rule {
            re: Regex::new(pattern)?,
            capture: 1,
            not_after,
        })
    }

    /// A user's regular expression, whose captures do not change its span.
    pub fn whole_match(pattern: &str) -> Result<Rule, regex::Error> {
        let mut rule = Self::new(pattern, None)?;
        rule.capture = 0;
        Ok(rule)
    }

    /// The spans of `re.finditer` over `line`, in order.
    pub fn find_iter<'a>(&'a self, line: &'a str) -> impl Iterator<Item = (usize, usize)> + 'a {
        let mut pos = 0;
        std::iter::from_fn(move || {
            while pos <= line.len() {
                let caps = self.re.captures_at(line, pos)?;
                let m = caps.get(self.capture).or_else(|| caps.get(0))?;
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

    /// `re.search(line) is not None`.
    pub fn is_match(&self, line: &str) -> bool {
        self.find_iter(line).next().is_some()
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
