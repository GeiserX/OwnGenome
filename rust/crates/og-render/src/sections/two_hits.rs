//! The genes a called deletion and a small variant both reach (og summary
//! collect's two-hit rule), shown under the carrier list: each gene with
//! its verdict, then which deletion callers were read.

use crate::model::Malformed;
use crate::sections::clinvar::D;

/// `str.capitalize()`: the first character upper case, the rest lower.
fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(first) => {
            pycompat::text::upper(&first.to_string()) + &pycompat::text::lower(c.as_str())
        }
        None => String::new(),
    }
}

/// `"; ".join(items)` of a list of strings.
fn joined(v: &D<'_>) -> Result<String, Malformed> {
    let items = v.items()?;
    let parts: Vec<&str> = items.iter().map(D::string).collect::<Result<_, _>>()?;
    Ok(parts.join("; "))
}

/// The lines, none for a summary written before the rule.
pub fn lines(t: Option<&D<'_>>) -> Result<Vec<String>, Malformed> {
    let Some(t) = t.filter(|t| t.truthy()) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::new();
    let genes = t.at("genes")?;
    for x in genes.items()? {
        let allele = match x.opt("allele")? {
            Some(a) => format!(" ({})", a.text()?),
            None => String::new(),
        };
        out.push(format!(
            "{}: {} ({}; AR), both copies may be hit: {} and {} {}{allele}; {}",
            capitalize(x.at("verdict")?.string()?),
            x.s("gene")?,
            x.s("condition")?,
            x.s("deletion")?,
            x.s("variant")?,
            x.s("genotype")?,
            x.s("why")?
        ));
    }
    let found = if genes.truthy() {
        format!("{} gene(s); ", genes.len()?)
    } else {
        "none found; ".to_owned()
    };
    let read = joined(&t.at("deletions_read")?)?;
    let read = if read.is_empty() {
        "no caller".to_owned()
    } else {
        read
    };
    let not_read = t.at("deletions_not_read")?;
    let not_read = if not_read.truthy() {
        format!("; not read: {}", joined(&not_read)?)
    } else {
        String::new()
    };
    out.push(format!(
        "Two hits across callers (a deletion plus a variant on the other copy): {found}deletions read from \
         {read}{not_read}"
    ));
    Ok(out)
}
