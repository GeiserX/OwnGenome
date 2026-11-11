use memchr::{memchr_iter, memmem};
use smallvec::SmallVec;

/// The fields of VEP's CSQ annotation, from its header line
/// `##INFO=<ID=CSQ,...,Description="... Format: Allele|Consequence|...">`.
///
/// Read as [`crate::rank`] reads a line that has one: the names after the last
/// `Format: `, up to the closing quote, split on `|`; a name given twice
/// means its last position. Indices start at 0.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CsqLayout {
    pub names: Vec<String>,
    pub impact: Option<usize>,
    pub mane: Option<usize>,
    pub canonical: Option<usize>,
    pub consequence: Option<usize>,
    pub symbol: Option<usize>,
    /// The frequency fields: `MAX_AF` when there is one, else `gnomADe_AF`
    /// and `gnomADg_AF`, those present.
    pub freq: SmallVec<[usize; 2]>,
    /// One past the highest index above: an entry needs splitting no further.
    pub max_needed: usize,
}

impl CsqLayout {
    /// Reads the CSQ header line (without its line ending); `None` for any
    /// other line, or a CSQ line with no `Format: `.
    pub fn from_header_line(line: &[u8]) -> Option<CsqLayout> {
        let rest = line.strip_prefix(b"##INFO=<ID=CSQ,")?;
        let at = memmem::rfind(rest, b"Format: ")?;
        let format = &rest[at + b"Format: ".len()..];
        let format = &format[..memchr::memchr(b'"', format).unwrap_or(format.len())];
        let names: Vec<String> = format
            .split(|&c| c == b'|')
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .collect();
        let index = |name: &str| names.iter().rposition(|n| n == name);
        let mut freq = SmallVec::new();
        if let Some(i) = index("MAX_AF") {
            freq.push(i);
        } else {
            freq.extend(index("gnomADe_AF"));
            freq.extend(index("gnomADg_AF"));
        }
        let mut layout = CsqLayout {
            impact: index("IMPACT"),
            mane: index("MANE_SELECT"),
            canonical: index("CANONICAL"),
            consequence: index("Consequence"),
            symbol: index("SYMBOL"),
            freq,
            max_needed: 0,
            names,
        };
        layout.max_needed = [
            layout.impact,
            layout.mane,
            layout.canonical,
            layout.consequence,
            layout.symbol,
        ]
        .into_iter()
        .flatten()
        .chain(layout.freq.iter().copied())
        .map(|i| i + 1)
        .max()
        .unwrap_or(0);
        Some(layout)
    }

    /// The position of field `name` (the last one when it is there twice).
    pub fn index(&self, name: &str) -> Option<usize> {
        self.names.iter().rposition(|n| n == name)
    }
}

/// Field `i` of one CSQ entry (split on `|`); empty when the entry has fewer
/// fields, as awk's `split` leaves a missing field empty.
pub fn csq_field(entry: &[u8], i: usize) -> &[u8] {
    let mut start = 0;
    let mut bars = memchr_iter(b'|', entry);
    for _ in 0..i {
        match bars.next() {
            Some(k) => start = k + 1,
            None => return &[],
        }
    }
    let end = bars.next().unwrap_or(entry.len());
    &entry[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout() {
        let line = br#"##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|CANONICAL|MANE_SELECT|gnomADe_AF|gnomADg_AF">"#;
        let l = CsqLayout::from_header_line(line).unwrap();
        assert_eq!(l.names.len(), 8);
        assert_eq!(
            (l.consequence, l.impact, l.symbol, l.canonical, l.mane),
            (Some(1), Some(2), Some(3), Some(4), Some(5))
        );
        assert_eq!(l.freq.as_slice(), [6, 7]);
        assert_eq!(l.max_needed, 8);
        assert_eq!(l.index("Allele"), Some(0));
        let max =
            br#"##INFO=<ID=CSQ,Description="x. Format: Allele|MAX_AF|IMPACT|MAX_AF|gnomADe_AF">"#;
        let l = CsqLayout::from_header_line(max).unwrap();
        assert_eq!(l.freq.as_slice(), [3], "MAX_AF wins, at its last position");
        assert_eq!((l.impact, l.mane, l.max_needed), (Some(2), None, 4));
        assert_eq!(
            CsqLayout::from_header_line(b"##INFO=<ID=CSQX,Description=\"Format: A\">"),
            None
        );
        assert_eq!(
            CsqLayout::from_header_line(b"##INFO=<ID=CSQ,Description=\"none\">"),
            None
        );
    }

    #[test]
    fn fields() {
        let e = b"A|missense_variant|MODERATE||YES";
        assert_eq!(csq_field(e, 0), b"A");
        assert_eq!(csq_field(e, 2), b"MODERATE");
        assert_eq!(csq_field(e, 3), b"");
        assert_eq!(csq_field(e, 4), b"YES");
        assert_eq!(csq_field(e, 5), b"");
        assert_eq!(csq_field(e, 50), b"");
        assert_eq!(csq_field(b"", 0), b"");
    }
}
