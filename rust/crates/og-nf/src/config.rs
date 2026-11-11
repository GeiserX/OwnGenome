//! `nextflow.config`: the plain literals of its `params { }` block, and its
//! `withName` selectors with their `errorStrategy` closures. Comments are
//! blanked first, strings respected.

use crate::groovy::{closing_brace, strip_comments};
use pycompat::json::Value;
use regex::Regex;
use std::sync::LazyLock;

static PARAMS: LazyLock<Regex> = LazyLock::new(|| re(r"(?ms)^params\s*\{(.*?)^\}"));
static PARAM: LazyLock<Regex> =
    LazyLock::new(|| re(r#"(?m)^\s*(\w+)\s*=\s*('[^']*'|"[^"]*"|[-\w.]+)"#));
static WITH_NAME: LazyLock<Regex> = LazyLock::new(|| re(r"withName:\s*'([^']*)'\s*\{"));
static STRATEGY: LazyLock<Regex> = LazyLock::new(|| re(r"\berrorStrategy\s*=\s*\{"));

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

/// `name = value` of the `params` block, in file order, for the values that
/// are plain literals: a quoted string, `true`, `false`, `null`, or what
/// Python's `int()` or else `float()` reads (so `1_000` and `nan` count). A
/// value it cannot read (an expression, a variable) is left out.
pub fn params(src: &str) -> Vec<(String, Value)> {
    let text = strip_comments(src);
    let Some(block) = PARAMS.captures(&text) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for c in PARAM.captures_iter(&block[1]) {
        let v = &c[2];
        let value = if v.starts_with(['\'', '"']) {
            Value::Str(v[1..v.len() - 1].to_owned())
        } else {
            match v {
                "true" => Value::Bool(true),
                "false" => Value::Bool(false),
                "null" => Value::Null,
                _ => match pycompat::num::parse_int(v) {
                    Ok(n) => Value::Int(n),
                    Err(_) => match pycompat::num::parse_float(v) {
                        Ok(x) => Value::Float(x),
                        Err(_) => continue,
                    },
                },
            }
        };
        out.push((c[1].to_owned(), value));
    }
    out
}

/// A `withName: '…' { … }` block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WithName {
    /// The selector, a regular expression over process names.
    pub selector: String,
    /// The body of its `errorStrategy = { … }` closure, when it has one.
    pub error_strategy: Option<String>,
}

/// Every `withName` block of the config, in file order. A block whose
/// braces never close runs to the end of the file.
pub fn with_names(src: &str) -> Vec<WithName> {
    let text = strip_comments(src);
    WITH_NAME
        .captures_iter(&text)
        .map(|c| {
            let m = c.get(0).expect("match");
            let open = m.end() - 1;
            let block = &text[m.end()..closing_brace(&text, open).unwrap_or(text.len())];
            let error_strategy = STRATEGY.find(block).map(|s| {
                let o = s.end() - 1;
                block[s.end()..closing_brace(block, o).unwrap_or(block.len())].to_owned()
            });
            WithName {
                selector: c[1].to_owned(),
                error_strategy,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_param_literals() {
        let src = "params {\n  a = 'x' // c\n  b = \"y//z\"\n  c = true\n  d = null\n  e = 1_000\n  f = 2.5\n  g = params.a\n  h = nan\n}\nother = 1\n";
        let p = params(src);
        let names: Vec<&str> = p.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, ["a", "b", "c", "d", "e", "f", "h"]);
        assert_eq!(p[1].1, Value::Str("y//z".into()));
        assert_eq!(p[2].1, Value::Bool(true));
        assert_eq!(p[3].1, Value::Null);
        assert_eq!(p[4].1, Value::Int(1000));
        assert_eq!(p[5].1, Value::Float(2.5));
        assert!(matches!(p[6].1, Value::Float(x) if x.is_nan()));
        assert!(params("nothing here").is_empty());
    }

    #[test]
    fn reads_with_name_blocks_with_nested_braces() {
        let src = "process {\n withName: 'A|B' { // 'ignore' here is a comment\n  cpus = 2\n  errorStrategy = { task.attempt == 1 ? 'retry' : { 'ignore' }() }\n }\n withName: 'C' { memory = '1 GB' }\n}\n";
        let w = with_names(src);
        assert_eq!(w.len(), 2);
        assert_eq!(w[0].selector, "A|B");
        let s = w[0].error_strategy.as_deref().unwrap();
        assert!(s.contains("'ignore'") && s.contains("'retry'"), "{s}");
        assert_eq!(w[1].selector, "C");
        assert_eq!(w[1].error_strategy, None);
    }
}
