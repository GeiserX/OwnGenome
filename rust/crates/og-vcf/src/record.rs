use std::fmt;

use memchr::{memchr, memchr_iter};

use crate::Gt;

/// Why a line is not a VCF record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// Fewer than the 8 fixed columns; holds the number found.
    Columns(usize),
    /// POS is not a whole number from 0 to 4294967295.
    Pos,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Columns(n) => write!(
                f,
                "a VCF record has 8 or more tab-separated columns, this line has {n}"
            ),
            Error::Pos => f.write_str("POS is not a whole number from 0 to 4294967295"),
        }
    }
}

impl std::error::Error for Error {}

/// The tabs a record keeps: the ends of CHROM to FORMAT.
const TABS: usize = 9;

/// One VCF record, read in place: every accessor returns a slice of the line.
#[derive(Debug, Clone, Copy)]
pub struct Record<'a> {
    line: &'a [u8],
    /// Offsets of the first tabs (at most [`TABS`]).
    tabs: [usize; TABS],
    /// How many of `tabs` were found.
    found: usize,
}

impl<'a> Record<'a> {
    /// Reads a line without its line ending. It needs the 8 fixed columns;
    /// nothing else is checked until it is asked for.
    pub fn parse(line: &'a [u8]) -> Result<Self, Error> {
        let mut tabs = [0; TABS];
        let mut found = 0;
        for (slot, at) in tabs.iter_mut().zip(memchr_iter(b'\t', line)) {
            *slot = at;
            found += 1;
        }
        if found < 7 {
            return Err(Error::Columns(found + 1));
        }
        Ok(Record { line, tabs, found })
    }

    /// The whole line.
    pub fn line(&self) -> &'a [u8] {
        self.line
    }

    /// Column `k` (0 is CHROM), for `k` up to 8 (FORMAT) when present.
    fn col(&self, k: usize) -> &'a [u8] {
        let start = if k == 0 { 0 } else { self.tabs[k - 1] + 1 };
        let end = if k < self.found {
            self.tabs[k]
        } else {
            self.line.len()
        };
        &self.line[start..end]
    }

    pub fn chrom(&self) -> &'a [u8] {
        self.col(0)
    }

    /// POS as written, digits only.
    pub fn pos(&self) -> Result<u32, Error> {
        let raw = self.col(1);
        if raw.is_empty() {
            return Err(Error::Pos);
        }
        raw.iter().try_fold(0u32, |n, &c| {
            if c.is_ascii_digit() {
                n.checked_mul(10)
                    .and_then(|n| n.checked_add(u32::from(c - b'0')))
                    .ok_or(Error::Pos)
            } else {
                Err(Error::Pos)
            }
        })
    }

    pub fn id(&self) -> &'a [u8] {
        self.col(2)
    }

    pub fn ref_(&self) -> &'a [u8] {
        self.col(3)
    }

    /// ALT as written (`.` when there is none).
    pub fn alt(&self) -> &'a [u8] {
        self.col(4)
    }

    /// The ALT alleles, split on `,`.
    pub fn alts(&self) -> impl Iterator<Item = &'a [u8]> + use<'a> {
        self.alt().split(|&c| c == b',')
    }

    pub fn qual(&self) -> &'a [u8] {
        self.col(5)
    }

    pub fn filter(&self) -> &'a [u8] {
        self.col(6)
    }

    /// INFO as written.
    pub fn info(&self) -> &'a [u8] {
        self.col(7)
    }

    /// INFO as `(key, value)` items in order; a flag has no value.
    pub fn info_items(&self) -> InfoItems<'a> {
        InfoItems::new(self.info())
    }

    /// The value of INFO `key`: `Some(None)` for a flag. When the key is
    /// there twice, the last one counts, as a Python dict or an awk array
    /// built from the items would keep it.
    pub fn info_get(&self, key: &[u8]) -> Option<Option<&'a [u8]>> {
        self.info_items()
            .filter(|(k, _)| *k == key)
            .last()
            .map(|(_, v)| v)
    }

    /// FORMAT, when the record has one.
    pub fn format(&self) -> Option<&'a [u8]> {
        (self.found >= 8).then(|| self.col(8))
    }

    /// The FORMAT keys in order (none without FORMAT).
    pub fn format_keys(&self) -> impl Iterator<Item = &'a [u8]> + use<'a> {
        self.format()
            .into_iter()
            .flat_map(|f| f.split(|&c| c == b':'))
    }

    /// The position of `key` in FORMAT.
    pub fn format_index(&self, key: &[u8]) -> Option<usize> {
        self.format_keys().position(|k| k == key)
    }

    /// Every sample column, as written.
    pub fn samples(&self) -> Option<&'a [u8]> {
        (self.found == TABS).then(|| &self.line[self.tabs[TABS - 1] + 1..])
    }

    /// Sample column `i` (0 is the first sample).
    pub fn sample(&self, i: usize) -> Option<&'a [u8]> {
        nth(self.samples()?, b'\t', i)
    }

    /// The value of FORMAT field `index` (from [`format_index`]) in sample
    /// `i`. `None` when the sample stops before it, as VCF allows.
    ///
    /// [`format_index`]: Record::format_index
    pub fn sample_value(&self, i: usize, index: usize) -> Option<&'a [u8]> {
        nth(self.sample(i)?, b':', index)
    }

    /// The genotype of sample `i`: GT must be the first FORMAT key.
    pub fn gt(&self, i: usize) -> Option<Gt> {
        if self.format_keys().next()? != b"GT" {
            return None;
        }
        Gt::parse(self.sample_value(i, 0)?)
    }

    /// Writes the line with `info` in place of INFO. With INFO rebuilt by
    /// [`write_info`] from [`info_items`](Record::info_items) unchanged, this
    /// writes the line as it was.
    pub fn write_with_info(&self, info: &[u8], out: &mut Vec<u8>) {
        out.extend_from_slice(&self.line[..self.tabs[6] + 1]);
        out.extend_from_slice(info);
        if self.found >= 8 {
            out.extend_from_slice(&self.line[self.tabs[7]..]);
        }
    }
}

/// Field `i` of `s` split on `sep`.
fn nth(s: &[u8], sep: u8, i: usize) -> Option<&[u8]> {
    let mut start = 0;
    let mut seps = memchr_iter(sep, s);
    for _ in 0..i {
        start = seps.next()? + 1;
    }
    let end = memchr(sep, &s[start..]).map_or(s.len(), |k| start + k);
    Some(&s[start..end])
}

/// INFO items in order. `.` alone has none; anything else, even an empty
/// field or an empty item between two `;`, gives back every item, so
/// [`write_info`] rebuilds the field byte for byte.
#[derive(Debug, Clone)]
pub struct InfoItems<'a> {
    rest: Option<&'a [u8]>,
}

impl<'a> InfoItems<'a> {
    pub fn new(info: &'a [u8]) -> Self {
        InfoItems {
            rest: (info != b".").then_some(info),
        }
    }
}

impl<'a> Iterator for InfoItems<'a> {
    type Item = (&'a [u8], Option<&'a [u8]>);

    fn next(&mut self) -> Option<Self::Item> {
        let s = self.rest?;
        let item = match memchr(b';', s) {
            Some(k) => {
                self.rest = Some(&s[k + 1..]);
                &s[..k]
            }
            None => {
                self.rest = None;
                s
            }
        };
        Some(match memchr(b'=', item) {
            Some(k) => (&item[..k], Some(&item[k + 1..])),
            None => (item, None),
        })
    }
}

/// Writes INFO items in the order given: `key=value` or a bare flag, joined
/// by `;`, and `.` when there are none.
pub fn write_info<'b>(
    items: impl IntoIterator<Item = (&'b [u8], Option<&'b [u8]>)>,
    out: &mut Vec<u8>,
) {
    let mut first = true;
    for (key, value) in items {
        if !first {
            out.push(b';');
        }
        first = false;
        out.extend_from_slice(key);
        if let Some(value) = value {
            out.push(b'=');
            out.extend_from_slice(value);
        }
    }
    if first {
        out.push(b'.');
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINE: &[u8] = b"chr2\t4294967295\trs1;rs2\tAC\tA,ACC\t12.5\tPASS;q10\tA=1;FLAG;B=x=y;;C=\tGT:AD:DP\t0/1:3,4:7\t1|2:1\t.";

    #[test]
    fn columns() {
        let r = Record::parse(LINE).unwrap();
        assert_eq!(r.chrom(), b"chr2");
        assert_eq!(r.pos(), Ok(u32::MAX));
        assert_eq!(r.id(), b"rs1;rs2");
        assert_eq!(r.ref_(), b"AC");
        assert_eq!(r.alts().collect::<Vec<_>>(), [&b"A"[..], b"ACC"]);
        assert_eq!(r.qual(), b"12.5");
        assert_eq!(r.filter(), b"PASS;q10");
        assert_eq!(r.format(), Some(&b"GT:AD:DP"[..]));
        assert_eq!(r.format_index(b"DP"), Some(2));
        assert_eq!(r.format_index(b"GQ"), None);
        assert_eq!(r.sample(0), Some(&b"0/1:3,4:7"[..]));
        assert_eq!(r.sample(1), Some(&b"1|2:1"[..]));
        assert_eq!(r.sample(2), Some(&b"."[..]));
        assert_eq!(r.sample(3), None);
        assert_eq!(r.sample_value(0, 1), Some(&b"3,4"[..]));
        assert_eq!(r.sample_value(1, 2), None, "a sample may stop early");
        assert_eq!(r.gt(1).unwrap().alleles.as_slice(), [Some(1), Some(2)]);
        assert!(r.gt(2).unwrap().alleles.as_slice() == [None]);
    }

    #[test]
    fn info_in_order_and_back() {
        let r = Record::parse(LINE).unwrap();
        let items: Vec<_> = r.info_items().collect();
        assert_eq!(
            items,
            [
                (&b"A"[..], Some(&b"1"[..])),
                (b"FLAG", None),
                (b"B", Some(b"x=y")),
                (b"", None),
                (b"C", Some(b"")),
            ]
        );
        assert_eq!(r.info_get(b"FLAG"), Some(None));
        assert_eq!(r.info_get(b"B"), Some(Some(&b"x=y"[..])));
        assert_eq!(r.info_get(b"Z"), None);
        let twice = Record::parse(b"1\t1\t.\tA\tC\t.\t.\tK=1;K=2").unwrap();
        assert_eq!(twice.info_get(b"K"), Some(Some(&b"2"[..])));
        for info in [&b"."[..], b"", b"A", b"A;", b";A", b"A=1;B", b";;"] {
            let mut out = Vec::new();
            write_info(InfoItems::new(info), &mut out);
            assert_eq!(out, info, "{}", String::from_utf8_lossy(info));
        }
    }

    #[test]
    fn rewrite_keeps_the_line() {
        for line in [
            LINE,
            b"1\t5\t.\tA\tC\t.\t.\t.",
            b"1\t5\t.\tA\tC\t.\t.\t",
            b"1\t5\t.\tA\tC\t.\t.\tX=1\tGT",
            b"1\t5\t.\tA\tC\t.\t.\tX=1\tGT\t0/0",
        ] {
            let r = Record::parse(line).unwrap();
            let mut info = Vec::new();
            write_info(r.info_items(), &mut info);
            let mut out = Vec::new();
            r.write_with_info(&info, &mut out);
            assert_eq!(out, line);
        }
        let r = Record::parse(b"1\t5\t.\tA\tC\t.\t.\tX=1").unwrap();
        assert_eq!(r.format(), None);
        assert_eq!(r.samples(), None);
        assert_eq!(r.gt(0), None);
    }

    #[test]
    fn errors() {
        assert_eq!(
            Record::parse(b"1\t5\t.\tA\tC\t.\t.").unwrap_err(),
            Error::Columns(7)
        );
        assert_eq!(Record::parse(b"").unwrap_err(), Error::Columns(1));
        for pos in ["", "-1", "+1", "1.0", " 1", "4294967296", "x"] {
            let line = format!("1\t{pos}\t.\tA\tC\t.\t.\t.");
            assert_eq!(
                Record::parse(line.as_bytes()).unwrap().pos(),
                Err(Error::Pos),
                "{pos:?}"
            );
        }
        assert_eq!(
            Record::parse(b"1\t0\t.\tA\tC\t.\t.\t.").unwrap().pos(),
            Ok(0)
        );
    }
}
