//! The GENCODE GTF reading the callability and SV rule builders share.
//!
//! A GTF line is read as bytes: its columns are found by their tabs, and an
//! attribute is found by a scan of the `key "value";` pairs of the ninth
//! column, only on the lines a builder keeps. GENCODE writes every pair as
//! `key "value"` or `key value` (a bare number), separated by `; `, and no
//! value holds a `;` or a `"`.

use super::BuildError;

/// One record of a GTF, by its columns.
pub(crate) struct Record<'a> {
    pub seq: &'a [u8],
    pub feature: &'a [u8],
    start: &'a [u8],
    end: &'a [u8],
    pub attrs: &'a [u8],
}

impl<'a> Record<'a> {
    /// The record of a data line (without its line ending); `None` for a
    /// comment.
    pub fn parse(line: &'a [u8]) -> Result<Option<Record<'a>>, BuildError> {
        if line.first() == Some(&b'#') {
            return Ok(None);
        }
        let mut c = line.splitn(9, |&b| b == b'\t');
        let mut next = || {
            c.next().ok_or_else(|| {
                BuildError::Input(format!(
                    "GTF: a line with under 9 columns: {:?}",
                    String::from_utf8_lossy(line)
                ))
            })
        };
        let seq = next()?;
        let _source = next()?;
        let feature = next()?;
        let start = next()?;
        let end = next()?;
        let (_score, _strand, _frame) = (next()?, next()?, next()?);
        let attrs = next()?;
        Ok(Some(Record {
            seq,
            feature,
            start,
            end,
            attrs,
        }))
    }

    /// The 0-based start and the end.
    pub fn span(&self) -> Result<(i64, i64), BuildError> {
        Ok((number(self.start)? - 1, number(self.end)?))
    }

    /// The value of attribute `key`; the last one when it is given twice,
    /// as a dict of the pairs keeps it. A value that is not UTF-8 is an
    /// error.
    pub fn attr(&self, key: &str) -> Result<Option<&'a str>, BuildError> {
        attr(self.attrs, key.as_bytes())
            .map(|v| {
                std::str::from_utf8(v)
                    .map_err(|_| BuildError::Input(format!("GTF: {key} is not UTF-8: {v:?}")))
            })
            .transpose()
    }
}

/// A position column.
fn number(cell: &[u8]) -> Result<i64, BuildError> {
    std::str::from_utf8(cell)
        .ok()
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| {
            BuildError::Input(format!(
                "GTF: not a whole number: {:?}",
                String::from_utf8_lossy(cell)
            ))
        })
}

/// The value of the last `key "value"` pair of an attribute column.
pub(crate) fn attr<'a>(attrs: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    let mut found = None;
    for pair in attrs.split(|&b| b == b';') {
        let pair = pair.trim_ascii_start();
        if let Some(rest) = pair.strip_prefix(key)
            && let Some(quoted) = rest.strip_prefix(b" \"")
            && let Some(end) = quoted.iter().position(|&b| b == b'"')
        {
            found = Some(&quoted[..end]);
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENE: &[u8] = b"chr1\tHAVANA\tgene\t11\t20\t.\t+\t.\tgene_id \"ENSG1\"; gene_type \"protein_coding\"; gene_name \"A B\"; level 2; tag \"x\"; tag \"y\"; xgene_name \"no\";";

    #[test]
    fn a_gene_line_reads_its_columns_and_attributes() {
        let r = Record::parse(GENE).unwrap().unwrap();
        assert_eq!((r.seq, r.feature), (&b"chr1"[..], &b"gene"[..]));
        assert_eq!(r.span().unwrap(), (10, 20));
        assert_eq!(r.attr("gene_name").unwrap(), Some("A B"), "not xgene_name");
        assert_eq!(
            r.attr("tag").unwrap(),
            Some("y"),
            "the last of a repeated key"
        );
        assert_eq!(
            r.attr("level").unwrap(),
            None,
            "an unquoted value is no pair"
        );
        assert_eq!(r.attr("gene").unwrap(), None, "a key is matched whole");
        assert_eq!(r.attr("hgnc_id").unwrap(), None);
    }

    #[test]
    fn comments_pass_and_short_lines_stop() {
        assert!(Record::parse(b"##description").unwrap().is_none());
        assert!(Record::parse(b"chr1\tx\tgene").is_err());
        let r = Record::parse(b"chr1\tx\tgene\tone\t2\t.\t+\t.\t")
            .unwrap()
            .unwrap();
        assert!(r.span().is_err());
    }
}
