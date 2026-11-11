//! The CHIP gene table, `assets/chip_genes.tsv`: which changes of each gene
//! count (`og chip check-table` checks it).
//!
//! Columns: `gene truncating missense note pmids`, tab separated; `#` lines
//! and blank lines are skipped. `truncating` is `all`, `none` or
//! `exons:A-B`; `missense` is `none` or a comma list of `all`, a codon
//! (`R882`), a change (`V617F`), a codon range (`102-292`) or `exons:A-B`.
//!
//! Numbers are ASCII digits only. The Python's `\d` also took other
//! Unicode digits; no table has them.

use crate::{Error, each_line};
use pycompat::num::py_strip;

/// The header the table must have.
pub const TABLE_COLUMNS: [&str; 5] = ["gene", "truncating", "missense", "note", "pmids"];

/// The 20 amino acids a missense rule may name.
const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// One gene of the table, as written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub gene: String,
    pub truncating: String,
    pub missense: String,
    pub note: String,
    pub pmids: String,
}

/// Which truncating changes count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Truncating {
    All,
    None,
    /// In exons A to B (ranks, inclusive).
    Exons(i128, i128),
}

/// One rule of the missense column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Missense {
    All,
    /// Any change of residue `ref` at codon `n` (`R882`).
    Codon(u8, i128),
    /// One change (`V617F`): residue, codon, new residue.
    Change(u8, i128, u8),
    /// Codons A to B.
    Range(i128, i128),
    /// Exons A to B.
    Exons(i128, i128),
}

/// A run of ASCII digits as a number (saturating: only order matters).
pub(crate) fn digits(s: &str) -> Option<i128> {
    if s.is_empty() || !s.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    Some(s.bytes().fold(0i128, |v, c| {
        v.saturating_mul(10).saturating_add(i128::from(c - b'0'))
    }))
}

/// `A-B` with 1 <= A <= B.
fn parse_range(tok: &str) -> Option<(i128, i128)> {
    let (a, b) = tok.split_once('-')?;
    let (a, b) = (digits(a)?, digits(b)?);
    (a <= b && a >= 1).then_some((a, b))
}

/// The truncating column; `None` when malformed.
pub fn parse_truncating(spec: &str) -> Option<Truncating> {
    match spec {
        "all" => Some(Truncating::All),
        "none" => Some(Truncating::None),
        _ => spec
            .strip_prefix("exons:")
            .and_then(parse_range)
            .map(|(a, b)| Truncating::Exons(a, b)),
    }
}

/// The missense column as rules (empty for `none`); `None` when a token is
/// malformed.
pub fn parse_missense(spec: &str) -> Option<Vec<Missense>> {
    if spec == "none" {
        return Some(Vec::new());
    }
    spec.split(',').map(parse_rule).collect()
}

fn parse_rule(tok: &str) -> Option<Missense> {
    if tok == "all" {
        return Some(Missense::All);
    }
    if let Some(r) = tok.strip_prefix("exons:") {
        return parse_range(r).map(|(a, b)| Missense::Exons(a, b));
    }
    if let Some((a, b)) = tok.split_once('-')
        && digits(a).is_some()
        && digits(b).is_some()
    {
        return parse_range(tok).map(|(a, b)| Missense::Range(a, b));
    }
    // [A-Z] digits [A-Z]?, each letter an amino acid.
    let b = tok.as_bytes();
    let first = *b.first()?;
    let last = *b.last()?;
    if !first.is_ascii_uppercase() || !AA.contains(&first) {
        return None;
    }
    let (body, alt) = if b.len() > 1 && last.is_ascii_uppercase() {
        (&tok[1..tok.len() - 1], Some(last))
    } else {
        (&tok[1..], None)
    };
    let n = digits(body)?;
    match alt {
        Some(a) if !AA.contains(&a) => None,
        Some(a) => Some(Missense::Change(first, n, a)),
        None => Some(Missense::Codon(first, n)),
    }
}

/// The rows of the table, or the problem with its layout.
pub fn read_table(path: &str) -> Result<Vec<Row>, Error> {
    let mut lines = Vec::new();
    each_line(path, |ln| {
        if !ln.starts_with('#') && !py_strip(ln).is_empty() {
            lines.push(ln.trim_end_matches('\n').to_owned());
        }
        Ok(())
    })?;
    let header = format!(
        "{path}: the header must be {} (tab separated)",
        TABLE_COLUMNS.join(" ")
    );
    match lines.first() {
        Some(h) if h.split('\t').eq(TABLE_COLUMNS) => {}
        _ => return Err(Error(header)),
    }
    lines[1..]
        .iter()
        .map(|ln| {
            let f: Vec<&str> = ln.split('\t').collect();
            match f[..] {
                [gene, truncating, missense, note, pmids] => Ok(Row {
                    gene: gene.to_owned(),
                    truncating: truncating.to_owned(),
                    missense: missense.to_owned(),
                    note: note.to_owned(),
                    pmids: pmids.to_owned(),
                }),
                _ => Err(Error(format!(
                    "{path}: '{}' has {} columns, not {}",
                    ln.chars().take(60).collect::<String>(),
                    f.len(),
                    TABLE_COLUMNS.len()
                ))),
            }
        })
        .collect()
}

fn is_symbol(g: &str) -> bool {
    !g.is_empty()
        && g.bytes()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == b'-')
}

fn is_pmids(s: &str) -> bool {
    s.split(';')
        .all(|p| (6..=9).contains(&p.len()) && p.bytes().all(|c| c.is_ascii_digit()))
}

/// The problems of the table, one line each (empty when it is fine).
pub fn check_table(path: &str) -> Result<Vec<String>, Error> {
    let mut problems = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for r in read_table(path)? {
        let g = &r.gene;
        if !is_symbol(g) {
            problems.push(format!("{g}: not a gene symbol"));
        }
        if !seen.insert(g.clone()) {
            problems.push(format!("{g}: listed twice"));
        }
        if parse_truncating(&r.truncating).is_none() {
            problems.push(format!(
                "{g}: truncating '{}' is not all, none or exons:A-B",
                r.truncating
            ));
        }
        if parse_missense(&r.missense).is_none() {
            problems.push(format!(
                "{g}: missense '{}' has a token that is not none, all, R882, V617F, 102-292 or exons:A-B",
                r.missense
            ));
        }
        if r.truncating == "none" && r.missense == "none" {
            problems.push(format!("{g}: counts no class at all"));
        }
        if py_strip(&r.note).is_empty() {
            problems.push(format!("{g}: no note"));
        }
        if !is_pmids(&r.pmids) {
            problems.push(format!(
                "{g}: pmids '{}' is not a ;-separated list of PubMed ids",
                r.pmids
            ));
        }
    }
    if seen.is_empty() {
        problems.push("the table lists no gene".to_owned());
    }
    Ok(problems)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truncating_specs() {
        assert_eq!(parse_truncating("all"), Some(Truncating::All));
        assert_eq!(parse_truncating("none"), Some(Truncating::None));
        assert_eq!(parse_truncating("exons:2-6"), Some(Truncating::Exons(2, 6)));
        for bad in ["exons:9-2", "exons:0-2", "exons:a-2", "exons:", "some", ""] {
            assert_eq!(parse_truncating(bad), None, "{bad}");
        }
    }

    #[test]
    fn missense_rules() {
        assert_eq!(
            parse_missense("R882,V617F,102-292,exons:5-8,all"),
            Some(vec![
                Missense::Codon(b'R', 882),
                Missense::Change(b'V', 617, b'F'),
                Missense::Range(102, 292),
                Missense::Exons(5, 8),
                Missense::All
            ])
        );
        assert_eq!(parse_missense("none"), Some(vec![]));
        for bad in [
            "R88X2", "B12", "R12B", "12", "292-102", "R", "", "R882,", "r882", "R882FF",
        ] {
            assert_eq!(parse_missense(bad), None, "{bad}");
        }
    }

    #[test]
    fn pmids_and_symbols() {
        assert!(is_pmids("25426837;15858187"));
        assert!(
            !is_pmids("PMID1") && !is_pmids("12345") && !is_pmids("1234567890") && !is_pmids("")
        );
        assert!(is_symbol("MLL-T3") && !is_symbol("dnmt3a") && !is_symbol(""));
    }
}
