//! The joint VCF as `og family report` and `og family segments` read it:
//! the columns of a record found once, each sample's genotype parsed at
//! most once per record.

use memchr::memchr_iter;
use og_pedigree::Die;
use smallvec::SmallVec;

/// Allele indices of a called genotype.
pub type Alleles = SmallVec<[u32; 2]>;

/// The tab-separated columns of a record.
pub type Cols<'a> = SmallVec<[&'a [u8]; 16]>;

/// Splits a record into its columns.
pub fn columns(line: &[u8]) -> Cols<'_> {
    let mut out = Cols::new();
    let mut start = 0;
    for t in memchr_iter(b'\t', line) {
        out.push(&line[start..t]);
        start = t + 1;
    }
    out.push(&line[start..]);
    out
}

/// A GT value as the Python read it: `(alleles, phased)`, alleles `None`
/// when the value is empty, `.`, or any allele is `.` or empty. `phased`
/// is whether a `|` appears. An allele that is not a whole number from 0
/// up is an error (the Python stopped with a traceback).
pub fn parse_gt(s: &[u8]) -> Result<(Option<Alleles>, bool), String> {
    if s.is_empty() || s == b"." {
        return Ok((None, false));
    }
    let phased = s.contains(&b'|');
    let parts = s.split(|&c| c == b'/' || c == b'|');
    let mut out = Alleles::new();
    for p in parts.clone() {
        if p.is_empty() || p == b"." {
            return Ok((None, phased));
        }
    }
    for p in parts {
        let v = std::str::from_utf8(p)
            .ok()
            .and_then(|t| pycompat::num::parse_int(t).ok())
            .and_then(|v| u32::try_from(v).ok())
            .ok_or_else(|| {
                format!(
                    "genotype '{}' is not a list of allele numbers",
                    String::from_utf8_lossy(s)
                )
            })?;
        out.push(v);
    }
    Ok((Some(out), phased))
}

/// `set(alleles)` as a short sorted list without repeats.
pub fn set(a: &[u32]) -> Alleles {
    let mut s: Alleles = a.iter().copied().collect();
    s.sort_unstable();
    s.dedup();
    s
}

/// The `(from_father, from_mother)` assignments a diploid child genotype
/// allows; a parent with no call (`None`) allows any allele.
pub fn resolve(child: &[u32], f: Option<&[u32]>, m: Option<&[u32]>) -> SmallVec<[(u32, u32); 2]> {
    let (a, b) = (child[0], child[1]);
    let mut opts = SmallVec::new();
    for (x, y) in [(a, b), (b, a)] {
        if f.is_none_or(|f| f.contains(&x))
            && m.is_none_or(|m| m.contains(&y))
            && !opts.contains(&(x, y))
        {
            opts.push((x, y));
        }
    }
    opts
}

/// One sample's column of the current record.
struct Entry<'a> {
    gt: Option<Alleles>,
    phased: bool,
    vals: SmallVec<[&'a [u8]; 8]>,
}

/// The samples of one record, each split and parsed on first use.
pub struct Calls<'a> {
    cols: &'a [&'a [u8]],
    gi: usize,
    cache: SmallVec<[Option<Entry<'a>>; 8]>,
    /// For messages: the file and line number.
    pub path: &'a str,
    pub line: u64,
}

impl<'a> Calls<'a> {
    /// `gi`: the index of GT in FORMAT.
    pub fn new(cols: &'a [&'a [u8]], gi: usize, path: &'a str, line: u64) -> Self {
        let mut cache = SmallVec::new();
        cache.resize_with(cols.len().saturating_sub(9), || None);
        Calls {
            cols,
            gi,
            cache,
            path,
            line,
        }
    }

    fn entry(&mut self, col: usize) -> Result<&Entry<'a>, Die> {
        let k = col - 9;
        if k >= self.cache.len() {
            return Err(Die::one(format!(
                "{} line {}: {} columns, the #CHROM line has more",
                self.path,
                self.line,
                self.cols.len()
            )));
        }
        if self.cache[k].is_none() {
            let vals: SmallVec<[&'a [u8]; 8]> = self.cols[col].split(|&c| c == b':').collect();
            let raw = vals.get(self.gi).copied().unwrap_or(b".");
            let (gt, phased) = parse_gt(raw)
                .map_err(|e| Die::one(format!("{} line {}: {e}", self.path, self.line)))?;
            self.cache[k] = Some(Entry { gt, phased, vals });
        }
        Ok(self.cache[k].as_ref().expect("filled above"))
    }

    /// The alleles of the sample in column `col`, `None` when not called.
    pub fn gt(&mut self, col: usize) -> Result<Option<Alleles>, Die> {
        Ok(self.entry(col)?.gt.clone())
    }

    /// The alleles and whether the GT is phased.
    pub fn gt_phased(&mut self, col: usize) -> Result<(Option<Alleles>, bool), Die> {
        let e = self.entry(col)?;
        Ok((e.gt.clone(), e.phased))
    }

    /// FORMAT value `i` of the sample in column `col`, `None` past its end.
    pub fn val(&mut self, col: usize, i: usize) -> Result<Option<&'a [u8]>, Die> {
        Ok(self.entry(col)?.vals.get(i).copied())
    }
}

/// The joint VCF is sorted, as the chromosome map needs: no position goes
/// back on a contig, and no contig comes back after another.
#[derive(Default)]
pub struct Sorted {
    prev: Option<(Vec<u8>, i64)>,
    done: std::collections::HashSet<Vec<u8>>,
}

impl Sorted {
    /// Checks the next record; the error names the record and the one before.
    pub fn next(&mut self, path: &str, chrom: &[u8], pos: i64) -> Result<(), Die> {
        let bad = match &self.prev {
            Some((pc, pp)) if pc.as_slice() == chrom => pos < *pp,
            Some((pc, _)) => {
                self.done.insert(pc.clone());
                self.done.contains(chrom)
            }
            None => false,
        };
        if bad {
            let (pc, pp) = self.prev.as_ref().expect("bad needs a record before");
            return Err(Die::one(format!(
                "{path} is not sorted: {}:{pos} after {}:{pp}",
                String::from_utf8_lossy(chrom),
                String::from_utf8_lossy(pc)
            )));
        }
        match &mut self.prev {
            Some((pc, pp)) if pc.as_slice() == chrom => *pp = pos,
            _ => self.prev = Some((chrom.to_vec(), pos)),
        }
        Ok(())
    }
}

/// POS as the Python's `int()` read it.
pub fn pos(raw: &[u8], path: &str, line: u64) -> Result<i64, Die> {
    std::str::from_utf8(raw)
        .ok()
        .and_then(|t| pycompat::num::parse_int(t).ok())
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| {
            Die::one(format!(
                "{path} line {line}: POS '{}' is not a whole number",
                String::from_utf8_lossy(raw)
            ))
        })
}

/// The index of `key` in a FORMAT column.
pub fn format_index(format: &[u8], key: &[u8]) -> Option<usize> {
    format.split(|&c| c == b':').position(|k| k == key)
}

/// The sample names of a `#CHROM` line, in column order.
pub fn sample_names(line: &[u8]) -> Vec<String> {
    columns(line)
        .iter()
        .skip(9)
        .map(|c| String::from_utf8_lossy(c).into_owned())
        .collect()
}

/// `{s: 9 + i}`: a name given twice maps to its last column.
pub fn sample_index(names: &[String]) -> std::collections::HashMap<String, usize> {
    names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), 9 + i))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_as_python_reads_it() {
        assert_eq!(
            parse_gt(b"0/1").unwrap(),
            (Some(Alleles::from_slice(&[0, 1])), false)
        );
        assert_eq!(
            parse_gt(b"1|0").unwrap(),
            (Some(Alleles::from_slice(&[1, 0])), true)
        );
        assert_eq!(parse_gt(b".|1").unwrap(), (None, true));
        assert_eq!(parse_gt(b".").unwrap(), (None, false));
        assert_eq!(parse_gt(b"").unwrap(), (None, false));
        assert_eq!(
            parse_gt(b"1").unwrap(),
            (Some(Alleles::from_slice(&[1])), false)
        );
        assert!(parse_gt(b"A/1").is_err());
    }

    #[test]
    fn resolve_as_python() {
        assert_eq!(
            resolve(&[0, 1], Some(&[0]), Some(&[0, 1])).as_slice(),
            &[(0, 1)]
        );
        assert_eq!(resolve(&[0, 1], Some(&[0, 1]), Some(&[0, 1])).len(), 2);
        assert_eq!(resolve(&[1, 1], None, Some(&[1])).as_slice(), &[(1, 1)]);
        assert!(resolve(&[1, 1], Some(&[0]), Some(&[0])).is_empty());
    }

    #[test]
    fn sorted_check() {
        let mut s = Sorted::default();
        assert!(s.next("v", b"chr1", 5).is_ok());
        assert!(s.next("v", b"chr1", 5).is_ok());
        assert!(s.next("v", b"chr2", 1).is_ok());
        let e = s.next("v", b"chr2", 0).unwrap_err();
        assert_eq!(e.lines, ["v is not sorted: chr2:0 after chr2:1"]);
        let mut s = Sorted::default();
        s.next("v", b"chr1", 5).unwrap();
        s.next("v", b"chr2", 1).unwrap();
        let e = s.next("v", b"chr1", 9).unwrap_err();
        assert_eq!(e.lines, ["v is not sorted: chr1:9 after chr2:1"]);
    }
}
