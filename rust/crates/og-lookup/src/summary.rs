//! What the reports (`summary.json`) say about a gene or an rsID: the same
//! names rule `bin/results_index.py` builds the search index with
//! (`names_of` in `bin/genome_lookup.py`), so a lookup and the search find
//! the same mentions.

use crate::Error;
use pycompat::json::{JsonError, Value};
use pycompat::text;
use pycompat::value::{py_str, truthy};
use std::collections::BTreeSet;

/// Keys of a summary entry that name a gene. "locus" is how the repeat
/// expansions, the CFTR/UGT1A1/SLC6A4 readings and the Locityper rows name
/// theirs (HTT, CFTR, HLA-DRB1).
const GENE_KEYS: [&str; 6] = ["gene", "genes", "Gene", "SYMBOL", "symbol", "locus"];

/// A loaded `summary.json`, or why it could not be read.
pub type Summary = Result<Value, Error>;

/// Reads `path` as Python's `json.load(open(path))` does (UTF-8, strict).
pub fn load(path: &str) -> Summary {
    let bytes = std::fs::read(path).map_err(|e| Error::io(&e, path))?;
    let text =
        String::from_utf8(bytes).map_err(|_| Error::value(format!("{path}: not UTF-8 text")))?;
    pycompat::json::loads(&text).map_err(|e| match e {
        JsonError::Decode(msg, at) => Error::value(decode_error(&text, msg, at)),
        other => Error::value(other.to_string()),
    })
}

/// A JSONDecodeError as Python words it: `MSG: line L column C (char N)`,
/// counted in characters.
fn decode_error(doc: &str, msg: &str, byte: usize) -> String {
    let before = &doc[..doc.floor_char_boundary(byte)];
    let pos = before.chars().count();
    let line = before.matches('\n').count() + 1;
    let col = match before.rfind('\n') {
        Some(nl) => before[nl + 1..].chars().count() + 1,
        None => pos + 1,
    };
    format!("{msg}: line {line} column {col} (char {pos})")
}

/// `^rs\d+$` in any case, `$` also before a final `\n` (Python's `re`).
fn is_rsid(s: &str) -> bool {
    let s = s.strip_suffix('\n').unwrap_or(s);
    let mut cs = s.chars();
    matches!((cs.next(), cs.next()), (Some('r' | 'R'), Some('s' | 'S')))
        && s.chars().count() > 2
        && s.chars().skip(2).all(text::is_decimal)
}

/// The genes (upper-cased) and rsIDs (lower-cased) a dict names itself: a
/// gene key's string is split on `, ; & / ( )` and spaces and keeps the
/// words with an upper-case letter ("VDBP (GC)" gives VDBP and GC, "SMN
/// exon 7-8 deletion" gives SMN); a list gives its strings, a dict its keys.
fn names_of(obj: &Value) -> (BTreeSet<String>, BTreeSet<String>) {
    let mut genes = BTreeSet::new();
    let mut rsids = BTreeSet::new();
    let Value::Object(o) = obj else {
        return (genes, rsids);
    };
    for k in GENE_KEYS {
        let names: Vec<&str> = match o.get(k) {
            Some(Value::Str(v)) => v
                .split([',', ';', '&', '/', ' ', '(', ')'])
                .filter(|n| *n != text::lower(n))
                .collect(),
            Some(Value::Array(a)) => a
                .iter()
                .filter_map(|x| match x {
                    Value::Str(s) => Some(s.as_str()),
                    _ => None,
                })
                .collect(),
            Some(Value::Object(d)) => d.keys().map(String::as_str).collect(),
            _ => continue,
        };
        genes.extend(
            names
                .into_iter()
                .filter(|n| !n.is_empty() && *n != ".")
                .map(text::upper),
        );
    }
    for k in ["rsid", "rsids"] {
        let vals: Vec<&Value> = match o.get(k) {
            Some(v @ Value::Str(_)) => vec![v],
            Some(Value::Array(a)) => a.iter().collect(),
            _ => Vec::new(),
        };
        for v in vals {
            if let Value::Str(s) = v
                && is_rsid(s)
            {
                rsids.insert(text::lower(s));
            }
        }
    }
    (genes, rsids)
}

/// Every dict in `obj` that names one of `genes` or `rsid`: the first dict
/// on each branch that names it, not the dicts inside it.
fn walk<'a>(
    obj: &'a Value,
    genes: &BTreeSet<String>,
    rsid: Option<&str>,
    out: &mut Vec<&'a Value>,
) {
    match obj {
        Value::Object(o) => {
            let (g, r) = names_of(obj);
            if !g.is_disjoint(genes) || rsid.is_some_and(|x| r.contains(x)) {
                out.push(obj);
                return;
            }
            for v in o.values() {
                walk(v, genes, rsid, out);
            }
        }
        Value::Array(a) => {
            for v in a {
                walk(v, genes, rsid, out);
            }
        }
        _ => {}
    }
}

/// A scalar Python's `isinstance(v, (str, int, float))` takes (bool is an int).
fn scalar(v: &Value) -> bool {
    matches!(
        v,
        Value::Str(_) | Value::Int(_) | Value::Float(_) | Value::Bool(_)
    )
}

/// `str(v)` of a scalar.
fn s(v: &Value) -> String {
    py_str(v).unwrap_or_default()
}

/// One entry as `key=value` pairs, cut at `limit` characters.
fn flat(d: &Value, limit: usize) -> String {
    let Value::Object(o) = d else {
        return String::new();
    };
    let mut parts = Vec::new();
    for (k, v) in o {
        match v {
            Value::Str(x) if x.is_empty() => {}
            v if scalar(v) => parts.push(format!("{k}={}", s(v))),
            Value::Array(a) if !a.is_empty() && a.iter().all(scalar) => {
                parts.push(format!(
                    "{k}={}",
                    a.iter().map(s).collect::<Vec<_>>().join(",")
                ));
            }
            Value::Object(m) if !m.is_empty() && m.values().all(scalar) => {
                let items: Vec<String> = m.iter().map(|(a, b)| format!("{a}: {}", s(b))).collect();
                parts.push(format!("{k}={}", items.join("; ")));
            }
            _ => {}
        }
    }
    let all = parts.join(", ");
    if all.chars().count() <= limit {
        all
    } else {
        let mut cut: String = all.chars().take(limit - 3).collect();
        cut.push_str("...");
        cut
    }
}

/// The "in the reports" lines for `genes` and `rsid`.
pub fn lines(
    summary: Option<&Summary>,
    genes: &BTreeSet<String>,
    rsid: Option<&str>,
) -> Result<Vec<String>, Error> {
    let Some(summary) = summary else {
        return Ok(vec![
            "  no summary.json for this sample: run scripts/generate-report.sh <sample> (on Nextflow output, \
             GENOME_DIR=<outdir>) to add what the reports say"
                .to_owned(),
        ]);
    };
    let summary = summary.as_ref().map_err(Clone::clone)?;
    let want: BTreeSet<String> = genes
        .iter()
        .filter(|g| !g.is_empty())
        .map(|g| text::upper(g))
        .collect();
    let mut out = Vec::new();
    let not_dict = || Error::other("summary.json: a section is not an object");
    let Value::Object(top) = summary else {
        return Err(not_dict());
    };
    let sections: Vec<(&String, &Value)> = match top.get("sections") {
        Some(v) if !truthy(v) => Vec::new(),
        None => Vec::new(),
        Some(Value::Object(m)) => m.iter().collect(),
        Some(_) => return Err(not_dict()),
    };
    for (key, sec) in sections {
        let Value::Object(sec_o) = sec else {
            return Err(not_dict());
        };
        let state = sec_o.get("state");
        let stale = matches!(state, Some(Value::Str(x)) if x == "stale");
        if !(stale || matches!(state, Some(Value::Str(x)) if x == "ok")) {
            continue;
        }
        let data = sec_o.get("data").filter(|v| truthy(v));
        let mut hits = Vec::new();
        if let Some(d) = data {
            walk(d, &want, rsid, &mut hits);
        }
        if hits.is_empty() {
            continue;
        }
        let title = match sec_o.get("title") {
            Some(v) if truthy(v) => s(v),
            _ => key.clone(),
        };
        let step = sec_o.get("step").map_or_else(|| "None".to_owned(), s);
        let stale = if stale {
            " (stale: from an earlier run)"
        } else {
            ""
        };
        out.push(format!("  {title} (step {step}){stale}:"));
        for h in hits.iter().take(20) {
            out.push(format!("    - {}", flat(h, 300)));
        }
        if hits.len() > 20 {
            let source = match sec_o.get("source") {
                Some(v) if truthy(v) => s(v),
                _ => "the section".to_owned(),
            };
            out.push(format!(
                "    - ... and {} more in {source}",
                hits.len() - 20
            ));
        }
    }
    if out.is_empty() {
        let mut names: BTreeSet<String> = want;
        if let Some(r) = rsid {
            names.insert(r.to_owned());
        }
        let names: Vec<String> = names.into_iter().collect();
        let what = if names.is_empty() {
            "it".to_owned()
        } else {
            names.join(", ")
        };
        out.push(format!("  none of the reports mentions {what}"));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(t: &str) -> Value {
        pycompat::json::loads(t).unwrap()
    }

    #[test]
    fn names_read_like_python() {
        let (g, r) = names_of(&v(
            r#"{"gene": "VDBP (GC)", "genes": ["a", "."], "locus": {"HTT": 1}, "rsids": ["RS12", "x", 3], "rsid": "rs5\n"}"#,
        ));
        assert_eq!(
            g.into_iter().collect::<Vec<_>>(),
            ["A", "GC", "HTT", "VDBP"]
        );
        assert_eq!(r.into_iter().collect::<Vec<_>>(), ["rs12", "rs5\n"]);
        let (g, _) = names_of(&v(r#"{"gene": "SMN exon 7-8 deletion"}"#));
        assert_eq!(g.into_iter().collect::<Vec<_>>(), ["SMN"]);
    }

    #[test]
    fn a_decode_error_reads_as_python_words_it() {
        assert_eq!(
            decode_error("{\n  \"é\": x", "Expecting value", 9),
            "Expecting value: line 2 column 7 (char 8)"
        );
    }

    #[test]
    fn flat_writes_python_str() {
        let d = v(
            r#"{"a": "", "b": true, "c": 1.0, "d": [1, "x", false], "e": {"k": null}, "f": {"k": 2}, "g": [], "h": null}"#,
        );
        assert_eq!(flat(&d, 300), "b=True, c=1.0, d=1,x,False, f=k: 2");
        assert_eq!(flat(&v(r#"{"a": "ééééé"}"#), 6), "a=é...");
    }
}
