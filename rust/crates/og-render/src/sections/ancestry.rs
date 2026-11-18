//! What the ancestry chapter's sections share: reading a row of the
//! summary the way the Python report did (`r["key"]`, `f"{v:<35}"`, a list
//! walked with `for`), the percentages of the fine-grained ancestry, and the
//! links where a haplogroup can be explored further (steps 12 and 37).
//!
//! Used by the fine-grained ancestry, haplogroup, Y chromosome, archaic,
//! ancient and polygenic score sections.

use pycompat::fmt::{fixed, ljust, rjust};
use pycompat::json::Value;

use crate::escape::e;
use crate::html::offsite_link;
use crate::model::{Line, Malformed, get, need, num, py, truthy};
use crate::notes::EXPLORE_NOTE;

/// Every explore address but the docs one carries the haplogroup in it.
pub const EXPLORE_DOCS: &str = "https://github.com/GeiserX/OwnGenome/blob/main/docs/";

/// The pages under [`EXPLORE_DOCS`] the collector links, one per haplogroup.
const EXPLORE_DOC_PAGES: [&str; 2] = [
    "37-y-haplogroup.md#explore-your-branch-further",
    "12-mito-haplogroup.md#explore-your-branch-further",
];

/// The mark of an address that leaves this computer with the haplogroup in it.
pub const OFFSITE_HAPLOGROUP: &str = "opens a website; the address carries the haplogroup";

/// `str(v)` of an optional value: `None` when it is absent.
pub fn py_or(v: Option<&Value>, default: &str) -> Result<String, Malformed> {
    v.map_or(Ok(default.to_owned()), py)
}

/// `d.get(key, default)` as text: the default only when the key is absent.
pub fn get_text(d: &Value, key: &str, default: &str, path: &str) -> Result<String, Malformed> {
    py_or(get(d, key, path)?, default)
}

/// `str(d[key])`.
pub fn need_text(d: &Value, key: &str, path: &str) -> Result<String, Malformed> {
    py(need(d, key, path)?)
}

/// The items `for x in v` walks: a list's items, a string's characters, a
/// dict's keys. `None` and numbers cannot be walked.
pub fn items(v: &Value, path: &str) -> Result<Vec<Value>, Malformed> {
    match v {
        Value::Array(a) => Ok(a.clone()),
        Value::Str(s) => Ok(s.chars().map(|c| Value::Str(c.to_string())).collect()),
        Value::Object(o) => Ok(o.keys().map(|k| Value::Str(k.clone())).collect()),
        _ => Err(Malformed(format!("{path} is not a list"))),
    }
}

/// `d.get(key, [])` walked with `for`: nothing when the key is absent.
pub fn list(d: &Value, key: &str, path: &str) -> Result<Vec<Value>, Malformed> {
    match get(d, key, path)? {
        None => Ok(Vec::new()),
        Some(v) => items(v, &format!("{path}.{key}")),
    }
}

/// `len(v)`.
pub fn py_len(v: &Value, path: &str) -> Result<usize, Malformed> {
    match v {
        Value::Array(a) => Ok(a.len()),
        Value::Str(s) => Ok(s.chars().count()),
        Value::Object(o) => Ok(o.len()),
        _ => Err(Malformed(format!("{path} has no length"))),
    }
}

/// Which side `f"{v:<N}"` (left) or `f"{v:>N}"` (right) pads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Align {
    Left,
    Right,
}

/// `f"{v:<N}"` or `f"{v:>N}"`: a string or a number padded to `width`
/// characters. A bool formats as the number it is, as in Python; `None`, a
/// list or a dict cannot take a width.
pub fn align(v: &Value, width: usize, side: Align, path: &str) -> Result<String, Malformed> {
    let text = match v {
        Value::Str(s) => s.clone(),
        Value::Int(i) => i.to_string(),
        Value::Float(_) => py(v)?,
        Value::Bool(b) => u8::from(*b).to_string(),
        _ => return Err(Malformed(format!("{path} cannot be padded"))),
    };
    Ok(match side {
        Align::Left => ljust(&text, width),
        Align::Right => rjust(&text, width),
    })
}

/// `f"{d[key]:<N}"` and `f"{d[key]:>N}"`.
pub fn need_align(
    d: &Value,
    key: &str,
    width: usize,
    side: Align,
    path: &str,
) -> Result<String, Malformed> {
    align(need(d, key, path)?, width, side, &format!("{path}.{key}"))
}

/// A number as Python compares it with `>=`; anything else is an error.
pub fn number(v: &Value, path: &str) -> Result<f64, Malformed> {
    num(Some(v)).ok_or_else(|| Malformed(format!("{path} is not a number")))
}

/// A fraction as a whole percentage (`95%`), `.` for none.
pub fn share(v: Option<&Value>, path: &str) -> Result<String, Malformed> {
    match v {
        None | Some(Value::Null) => Ok(".".to_owned()),
        Some(x) => Ok(format!("{}%", fixed(100.0 * number(x, path)?, 0))),
    }
}

/// A stat line whose value is markup built here.
pub fn stat_raw(label: &str, value: &str) -> Line {
    Line::Raw(format!(
        "    <div class=\"stat\"><span class=\"label\">{}</span><span class=\"value\">{value}</span></div>",
        e(label)
    ))
}

/// One explore entry of a haplogroup section: its name and address.
fn explore(d: &Value, path: &str) -> Result<Vec<(String, Value)>, Malformed> {
    let at = format!("{path}.explore[]");
    list(d, "explore", path)?
        .iter()
        .map(|x| Ok((need_text(x, "name", &at)?, need(x, "url", &at)?.clone())))
        .collect()
}

/// The text report's explore lines: `name: url` for each link, then the
/// note that nothing is uploaded; nothing when there is no link.
pub fn explore_lines(d: &Value, path: &str) -> Result<Vec<String>, Malformed> {
    let mut out = Vec::new();
    for (name, url) in explore(d, path)? {
        out.push(format!("{name}: {}", py(&url)?));
    }
    if !out.is_empty() {
        out.push(EXPLORE_NOTE.to_owned());
    }
    Ok(out)
}

/// The HTML report's explore lines: one stat per link, each marked as
/// leaving this computer, and whether the card shows the explore caveat.
/// The address is escaped, so a branch name cannot open a tag or end the
/// attribute.
pub fn explore_html(d: &Value, path: &str) -> Result<(Vec<Line>, bool), Malformed> {
    let mut out = Vec::new();
    for (name, url) in explore(d, path)? {
        let Value::Str(url) = url else {
            return Err(Malformed(format!("{path}.explore[].url is not a text")));
        };
        // Only the two fixed pages the collector links get the plain
        // warning: any other address, even under the docs folder, may carry
        // the haplogroup in its query.
        let why = if url
            .strip_prefix(EXPLORE_DOCS)
            .is_some_and(|page| EXPLORE_DOC_PAGES.contains(&page))
        {
            crate::html::OFFSITE
        } else {
            OFFSITE_HAPLOGROUP
        };
        // Only a web address with a host becomes a link: a `javascript:` or
        // `data:` address in a summary would run when clicked, and a bare
        // `https://` leads nowhere, so both show as text.
        let web = ["http://", "https://"].iter().any(|s| {
            url.get(..s.len())
                .is_some_and(|head| head.eq_ignore_ascii_case(s))
                && url[s.len()..]
                    .split(['/', '?', '#'])
                    .next()
                    .is_some_and(|host| !host.is_empty())
        });
        let value = if web {
            offsite_link(&url, &url, why)
        } else {
            e(&url).into_owned()
        };
        out.push(stat_raw(&name, &value));
    }
    let note = !out.is_empty();
    Ok((out, note))
}

/// `d.get(key, ".") not in (".", "")`: a value worth a line.
pub fn shown<'a>(d: &'a Value, key: &str, path: &str) -> Result<Option<&'a Value>, Malformed> {
    Ok(get(d, key, path)?.filter(|v| !matches!(v, Value::Str(s) if s.is_empty() || s == ".")))
}

/// `bool(d.get(key))`.
pub fn has(d: &Value, key: &str, path: &str) -> Result<bool, Malformed> {
    Ok(get(d, key, path)?.is_some_and(truthy))
}
