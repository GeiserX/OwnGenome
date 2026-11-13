//! The records of one region of a VCF, through its tabix index when it has
//! one, else by streaming the whole file and keeping only that region.

use crate::Error;
use crate::query::norm_chrom;
use og_io::{Lines, Utf8, tabix};
use std::io::Write;
use std::path::Path;

/// Alleles that stand for "any other allele", not a variant.
const SYMBOLIC: [&str; 3] = ["<*>", "<NON_REF>", "*"];

/// The CSQ field names of VEP's header line, in order.
#[derive(Debug, Clone, Default)]
pub struct Csq(pub Vec<String>);

impl Csq {
    /// The names from the `##INFO=<ID=CSQ,...Format: a|b|c">` line among
    /// `header`, or none.
    pub fn from_header<'a>(header: impl IntoIterator<Item = &'a str>) -> Csq {
        for h in header {
            if let Some(rest) = h.strip_prefix("##INFO=<ID=CSQ,")
                && let Some(at) = rest.find("Format: ")
            {
                let tail = &rest[at + "Format: ".len()..];
                // `[^"]+"`: at least one character, then a quote.
                if let Some(q) = tail.find('"').filter(|&q| q > 0) {
                    let names = pycompat::text::strip(&tail[..q], None);
                    return Csq(names.split('|').map(str::to_owned).collect());
                }
            }
        }
        Csq(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// One CSQ entry: its `|`-separated values, read by field name as
/// `dict(zip(fields, values))` reads them (the last field of a repeated
/// name wins; a field past the end of the values is missing).
#[derive(Debug, Clone)]
pub struct Entry<'a> {
    fields: &'a Csq,
    values: Vec<&'a str>,
}

impl<'a> Entry<'a> {
    /// The value of `name`, or `""` when the entry has none.
    pub fn get(&self, name: &str) -> &'a str {
        let n = self.fields.0.len().min(self.values.len());
        (0..n)
            .rev()
            .find(|&i| self.fields.0[i] == name)
            .map_or("", |i| self.values[i])
    }
}

/// A data line of a VCF.
#[derive(Debug, Clone)]
pub struct Rec {
    /// The chromosome, normalised ([`norm_chrom`]).
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_: String,
    pub alts: Vec<String>,
    pub qual: String,
    pub filter: String,
    /// INFO items in order; a flag has no value.
    info: Vec<(String, Option<String>)>,
    /// FORMAT keys paired with the first sample's values.
    fmt: Vec<(String, String)>,
    /// The last base: INFO/END when it is a number, else the last base of REF.
    pub end: u64,
}

impl Rec {
    /// Reads the tab-separated columns of a line (8 at least).
    pub fn parse(c: &[&str]) -> Result<Rec, Error> {
        if c.len() < 8 {
            return Err(Error::other(format!(
                "a VCF line has {} columns, fewer than 8",
                c.len()
            )));
        }
        let pos = int(c[1])?;
        let pos = u64::try_from(pos)
            .map_err(|_| Error::value(format!("a negative position: {}", c[1])))?;
        let info: Vec<(String, Option<String>)> = c[7]
            .split(';')
            .map(|kv| match kv.split_once('=') {
                Some((k, v)) => (k.to_owned(), Some(v.to_owned())),
                None => (kv.to_owned(), None),
            })
            .collect();
        let fmt = if c.len() >= 10 {
            c[8].split(':')
                .zip(c[9].split(':'))
                .map(|(k, v)| (k.to_owned(), v.to_owned()))
                .collect()
        } else {
            Vec::new()
        };
        let mut r = Rec {
            chrom: norm_chrom(c[0]),
            pos,
            id: c[2].to_owned(),
            ref_: c[3].to_owned(),
            alts: if c[4] == "." {
                Vec::new()
            } else {
                c[4].split(',').map(str::to_owned).collect()
            },
            qual: c[5].to_owned(),
            filter: c[6].to_owned(),
            info,
            fmt,
            end: 0,
        };
        r.end = match r.info_str("END") {
            Some(e) if !e.is_empty() && e.bytes().all(|b| b.is_ascii_digit()) => {
                e.parse().unwrap_or(u64::MAX)
            }
            _ => pos + (r.ref_.chars().count().max(1) as u64) - 1,
        };
        Ok(r)
    }

    /// INFO `key`, when it has a value (a flag has none). The last item of
    /// a repeated key wins.
    pub fn info_str(&self, key: &str) -> Option<&str> {
        match self.info.iter().rev().find(|(k, _)| k == key) {
            Some((_, Some(v))) => Some(v),
            _ => None,
        }
    }

    /// The first sample's value of FORMAT `key`.
    pub fn fmt(&self, key: &str) -> Option<&str> {
        self.fmt
            .iter()
            .rev()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// The alternate alleles that are variants.
    pub fn real_alts(&self) -> Vec<&str> {
        self.alts
            .iter()
            .map(String::as_str)
            .filter(|a| !SYMBOLIC.contains(a))
            .collect()
    }

    /// A reference block: no alternate allele but symbolic ones.
    pub fn is_block(&self) -> bool {
        self.real_alts().is_empty()
    }

    /// GT's allele indexes (`None` for `.`), and whether it is phased.
    pub fn gt(&self) -> Result<(Vec<Option<i128>>, bool), Error> {
        let g = self.fmt("GT").unwrap_or("");
        if g.is_empty() {
            return Ok((Vec::new(), false));
        }
        let idx = g
            .split(['/', '|'])
            .map(|a| if a == "." { Ok(None) } else { int(a).map(Some) })
            .collect::<Result<_, _>>()?;
        Ok((idx, g.contains('|')))
    }

    /// GT as alleles (`None` for a no-call or an index past the alleles),
    /// and whether it is phased.
    pub fn bases(&self) -> Result<(Vec<Option<&str>>, bool), Error> {
        let (idx, phased) = self.gt()?;
        let n = 1 + self.alts.len() as i128;
        let allele = |i: i128| -> Result<Option<&str>, Error> {
            // Python's list index: a negative one counts from the end.
            let j = if i < 0 { n + i } else { i };
            if i >= n {
                Ok(None)
            } else if j < 0 {
                Err(Error::other("list index out of range"))
            } else if j == 0 {
                Ok(Some(self.ref_.as_str()))
            } else {
                Ok(Some(self.alts[(j - 1) as usize].as_str()))
            }
        };
        let bases = idx
            .into_iter()
            .map(|i| i.map_or(Ok(None), allele))
            .collect::<Result<_, _>>()?;
        Ok((bases, phased))
    }

    /// Whether GT calls a variant allele.
    pub fn has_alt_call(&self) -> Result<bool, Error> {
        let (idx, _) = self.gt()?;
        Ok(idx.into_iter().flatten().any(|i| {
            i > 0
                && i <= self.alts.len() as i128
                && !SYMBOLIC.contains(&self.alts[(i - 1) as usize].as_str())
        }))
    }

    /// The CSQ entries, read by `fields`.
    pub fn csq<'a>(&'a self, fields: &'a Csq) -> Vec<Entry<'a>> {
        match self.info_str("CSQ") {
            Some(raw) if !fields.is_empty() => raw
                .split(',')
                .map(|e| Entry {
                    fields,
                    values: e.split('|').collect(),
                })
                .collect(),
            _ => Vec::new(),
        }
    }

    /// Every ID of the record: its ID column and VEP's Existing_variation,
    /// sorted and once each.
    pub fn ids(&self, fields: &Csq) -> Vec<String> {
        let mut out: Vec<String> = self
            .id
            .split(';')
            .filter(|i| *i != ".")
            .map(str::to_owned)
            .collect();
        for e in self.csq(fields) {
            out.extend(
                e.get("Existing_variation")
                    .split('&')
                    .filter(|x| !x.is_empty())
                    .map(str::to_owned),
            );
        }
        out.sort();
        out.dedup();
        out
    }

    /// The SYMBOLs of the CSQ entries.
    pub fn symbols(&self, fields: &Csq) -> Vec<String> {
        let mut out: Vec<String> = self
            .csq(fields)
            .iter()
            .map(|e| e.get("SYMBOL"))
            .filter(|s| !s.is_empty())
            .map(str::to_owned)
            .collect();
        out.sort();
        out.dedup();
        out
    }
}

/// Python's `int(s)`, or the `ValueError` it raises.
pub fn int(s: &str) -> Result<i128, Error> {
    pycompat::num::parse_int(s).map_err(|_| {
        Error::value(format!(
            "invalid literal for int() with base 10: {}",
            crate::repr(s)
        ))
    })
}

/// The records of one region, with the header's CSQ fields.
#[derive(Debug, Default)]
pub struct Region {
    pub csq: Csq,
    pub records: Vec<Rec>,
}

/// A region to read: a normalised chromosome and 1-based inclusive bounds.
#[derive(Debug, Clone, Copy)]
pub struct Span<'a> {
    pub chrom: &'a str,
    pub start: u64,
    pub end: u64,
}

/// The records of `path` that overlap `span` (a reference block by its
/// END, a deletion by its REF), in file order: through the `.tbi` or `.csi`
/// beside the file when there is one, else by streaming the whole file,
/// which `note` is told.
pub fn read(path: &str, span: Span, note: &mut dyn Write) -> Result<Region, Error> {
    {
        let has_index = [".tbi", ".csi"]
            .iter()
            .any(|x| Path::new(&format!("{path}{x}")).is_file());
        if !has_index {
            let _ = writeln!(
                note,
                "NOTE: {path} has no .tbi index; reading the whole file (bcftools index -t makes this fast)."
            );
        } else {
            match read_indexed(path, span) {
                Ok(r) => return Ok(r),
                Err(e) => {
                    let _ = writeln!(
                        note,
                        "NOTE: the index of {path} could not be read for {}:{}-{} ({}); reading the whole file.",
                        crate::query::contig_name(span.chrom),
                        span.start,
                        span.end,
                        e.msg
                    );
                }
            }
        }
    }
    read_stream(path, span)
}

/// Keeps `r` when it overlaps `span`.
fn overlaps(r: &Rec, span: Span) -> bool {
    r.chrom == span.chrom && r.pos <= span.end && r.end >= span.start
}

/// The header through the top of the file, then each index sequence that
/// names the region's chromosome, through the index.
fn read_indexed(path: &str, span: Span) -> Result<Region, Error> {
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut header = Vec::new();
    for line in Lines::new(file, Utf8::Replace) {
        let line = line.map_err(|e| Error::io(&e, path))?;
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("##") {
            header.push(line);
        }
    }
    let csq = Csq::from_header(header.iter().map(String::as_str));
    let mut reader = tabix::Reader::open(path).map_err(|e| Error::io(&e, path))?;
    let names: Vec<String> = reader
        .index()
        .names()
        .iter()
        .filter(|n| norm_chrom(n) == span.chrom)
        .cloned()
        .collect();
    let mut records = Vec::new();
    for name in names {
        for line in reader.query(&name, span.start.saturating_sub(1), span.end) {
            let line = line.map_err(|e| Error::io(&e, path))?;
            let text = String::from_utf8_lossy(&line);
            // A lone \r ends a line in Python's text mode.
            for part in text.split('\r') {
                if let Some(r) = region_line(part, span)? {
                    records.push(r);
                }
            }
        }
    }
    Ok(Region { csq, records })
}

/// A data line of the region, or `None` for one outside it (or too short
/// to be a record).
fn region_line(line: &str, span: Span) -> Result<Option<Rec>, Error> {
    if line.starts_with('#') {
        return Ok(None);
    }
    let c: Vec<&str> = line.split('\t').collect();
    if c.len() < 8 {
        return Ok(None);
    }
    if norm_chrom(c[0]) != span.chrom
        || c[1].is_empty()
        || !c[1].bytes().all(|b| b.is_ascii_digit())
        || c[1].parse::<u64>().map_or(true, |p| p > span.end)
    {
        return Ok(None);
    }
    let r = Rec::parse(&c)?;
    Ok(overlaps(&r, span).then_some(r))
}

/// Streams the whole file (plain, gzip or BGZF on several cores).
fn read_stream(path: &str, span: Span) -> Result<Region, Error> {
    let file = og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(|e| Error::io(&e, path))?;
    let mut lines = Lines::new(file, Utf8::Replace);
    let mut header = Vec::new();
    let mut records = Vec::new();
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        let l = line.strip_suffix('\n').unwrap_or(&line);
        if l.starts_with("##") {
            header.push(l.to_owned());
        } else if !l.starts_with('#')
            && let Some(r) = region_line(l, span)?
        {
            records.push(r);
        }
    }
    let csq = Csq::from_header(header.iter().map(String::as_str));
    Ok(Region { csq, records })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(line: &str) -> Rec {
        Rec::parse(&line.split('\t').collect::<Vec<_>>()).unwrap()
    }

    #[test]
    fn a_record_ends_at_its_end_or_its_ref() {
        assert_eq!(rec("chr7\t100\t.\tA\t<*>\t0\t.\tEND=200\tGT\t0/0").end, 200);
        assert_eq!(rec("chr7\t500\t.\tACGT\tA\t0\t.\t.\tGT\t0/1").end, 503);
        assert_eq!(rec("chr7\t500\t.\tA\tG\t0\t.\tEND=x\tGT\t0/1").end, 500);
    }

    #[test]
    fn genotypes_read_as_python_read_them() {
        let r = rec("chr7\t1\t.\tC\tT,<*>\t.\t.\t.\tGT:DP\t1|-1:9");
        let (b, phased) = r.bases().unwrap();
        assert_eq!(b, vec![Some("T"), Some("<*>")]);
        assert!(phased);
        assert!(r.has_alt_call().unwrap());
        let r = rec("chr7\t1\t.\tC\tT\t.\t.\t.\tGT\t0/");
        assert_eq!(r.gt().unwrap_err().code, 2);
        let r = rec("chr7\t1\t.\tC\tT\t.\t.\t.\tGT\t0/5");
        assert_eq!(r.bases().unwrap().0, vec![Some("C"), None]);
    }

    #[test]
    fn csq_fields_read_like_a_dict_of_a_zip() {
        let f = Csq(vec!["A".into(), "B".into(), "A".into()]);
        let r = rec("chr7\t1\t.\tC\tT\t.\t.\tCSQ=x|y|z,p|q\tGT\t0/1");
        let e = r.csq(&f);
        assert_eq!((e[0].get("A"), e[0].get("B")), ("z", "y"));
        assert_eq!((e[1].get("A"), e[1].get("C")), ("p", ""));
        let h = Csq::from_header([
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP. Format: Allele|SYMBOL\">",
        ]);
        assert_eq!(h.0, ["Allele", "SYMBOL"]);
    }
}
