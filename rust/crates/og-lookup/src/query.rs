//! What a query asks for, and how contig names compare across files.

use crate::Error;
use pycompat::text;

/// A parsed query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Query {
    /// An rsID, lower-cased (`rs4988235`).
    Rsid(String),
    /// A GRCh38 position: the normalised chromosome and the 1-based position.
    Site(String, u64),
    /// A gene symbol as given.
    Gene(String),
}

/// `rs` and one or more ASCII digits, in any case.
pub fn is_rsid(q: &str) -> bool {
    let b = q.as_bytes();
    b.len() > 2 && b[..2].eq_ignore_ascii_case(b"rs") && b[2..].iter().all(u8::is_ascii_digit)
}

/// Reads a query: an rsID (`rs123`), a position (`chr7:117,559,590`,
/// `7_117559590`, `chrMT-5`) or a gene symbol (`CFTR`). Digits are ASCII
/// only; a comma anywhere after the first digit is dropped, so `chr7:1,2,3`
/// is position 123, as it always was.
pub fn parse(q: &str) -> Result<Query, Error> {
    let q = text::strip(q, None);
    if is_rsid(q) {
        return Ok(Query::Rsid(q.to_ascii_lowercase()));
    }
    if let Some(at) = q.find([':', '_', '-']) {
        let (left, digits) = (&q[..at], &q[at + 1..]);
        let chrom = match left.get(..3) {
            Some(p) if p.eq_ignore_ascii_case("chr") => &left[3..],
            _ => left,
        };
        let chrom_ok = matches!(chrom.len(), 1 | 2) && chrom.bytes().all(|c| c.is_ascii_digit())
            || ["X", "Y", "M", "MT"]
                .iter()
                .any(|c| c.eq_ignore_ascii_case(chrom));
        let pos_ok = digits.as_bytes().first().is_some_and(u8::is_ascii_digit)
            && digits.bytes().all(|c| c.is_ascii_digit() || c == b',');
        if chrom_ok && pos_ok {
            let clean: String = digits.chars().filter(|&c| c != ',').collect();
            let pos: u64 = clean
                .parse()
                .map_err(|_| Error::value(format!("position is too large: {q}")))?;
            if pos < 1 {
                return Err(Error::value(format!("position must be 1 or more: {q}")));
            }
            return Ok(Query::Site(norm_chrom(chrom), pos));
        }
    }
    let b = q.as_bytes();
    if b.first().is_some_and(u8::is_ascii_alphanumeric)
        && b.iter()
            .all(|&c| c.is_ascii_alphanumeric() || matches!(c, b'.' | b'_' | b'-'))
    {
        return Ok(Query::Gene(q.to_owned()));
    }
    Err(Error::value(format!(
        "not an rsID (rs123), a position (chr7:117559590) or a gene symbol (CFTR): {}",
        crate::repr(q)
    )))
}

/// `chr7` and `7` read as `7`; `chrM`, `M` and `MT` as `MT`; the rest
/// upper-cased.
pub fn norm_chrom(c: &str) -> String {
    let c = match c.get(..3) {
        Some(p) if text::lower(p) == "chr" => &c[3..],
        _ => c,
    };
    let up = text::upper(c);
    if up == "M" || up == "MT" {
        "MT".to_owned()
    } else {
        up
    }
}

/// The pipeline's GRCh38 contig name of a normalised chromosome.
pub fn contig_name(c: &str) -> String {
    if c == "MT" {
        "chrM".to_owned()
    } else {
        format!("chr{c}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_each_kind_of_query() {
        assert_eq!(parse("RS4988235").unwrap(), Query::Rsid("rs4988235".into()));
        assert_eq!(
            parse(" chr7:1,234 ").unwrap(),
            Query::Site("7".into(), 1234)
        );
        assert_eq!(parse("7:1234").unwrap(), Query::Site("7".into(), 1234));
        assert_eq!(parse("chrMT:5").unwrap(), Query::Site("MT".into(), 5));
        assert_eq!(parse("CHRx_5").unwrap(), Query::Site("X".into(), 5));
        assert_eq!(parse("chr7:1,2,3").unwrap(), Query::Site("7".into(), 123));
        assert_eq!(parse("CFTR").unwrap(), Query::Gene("CFTR".into()));
        assert_eq!(parse("HLA-A").unwrap(), Query::Gene("HLA-A".into()));
        assert_eq!(parse("chr7").unwrap(), Query::Gene("chr7".into()));
        assert_eq!(contig_name("MT"), "chrM");
    }

    #[test]
    fn refuses_what_is_none_of_them() {
        for bad in ["chr7:0", "a b", "rs12;x", "", "chr7:,1", "-x"] {
            assert!(parse(bad).is_err(), "{bad:?} was accepted");
        }
        // Unicode digits are not digits here (Python's \d took them).
        assert!(parse("rs\u{661}\u{662}").is_err());
        let e = parse("a b").unwrap_err();
        assert_eq!(e.code, 2);
        assert!(e.msg.ends_with("gene symbol (CFTR): 'a b'"), "{}", e.msg);
    }

    #[test]
    fn contigs_compare_across_namings() {
        assert_eq!(norm_chrom("chr7"), "7");
        assert_eq!(norm_chrom("CHRm"), "MT");
        assert_eq!(norm_chrom("x"), "X");
        assert_eq!(norm_chrom("chrUn_KI270302v1"), "UN_KI270302V1");
    }
}
