//! The tables `og cn` reads: `samtools bedcov` sums, `samtools mpileup`
//! bases columns and a BAM's header.

use crate::catalog::Locus;
use crate::engine::{Means, SiteCounts};
use crate::{Error, each_line, py_float, py_int, strip};
use og_filters::py_repr;
use std::collections::{HashMap, HashSet};

/// One `samtools bedcov` line: chrom, start, end, name, depth sum.
#[derive(Debug, Clone, PartialEq)]
pub struct BedcovLine {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub sum: f64,
}

/// Every line of a `samtools bedcov` table (blank and `#` lines skipped).
pub fn read_bedcov_all(path: &str) -> Result<Vec<BedcovLine>, Error> {
    let mut out = Vec::new();
    each_line(path, |line| {
        if strip(line).is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            return Err(Error(format!(
                "{path}: expected chrom, start, end, name, depth sum; got {}",
                py_repr(strip(line))
            )));
        }
        out.push(BedcovLine {
            chrom: c[0].to_owned(),
            start: py_int(c[1])?,
            end: py_int(c[2])?,
            name: c[3].to_owned(),
            sum: py_float(c[4])?,
        });
        Ok(())
    })?;
    Ok(out)
}

/// `{name: mean depth}` from `samtools bedcov` over `row`'s segments, which
/// must all be there and nothing else; empty for an empty table (no reads
/// on that contig).
pub fn read_bedcov(path: &str, row: &Locus) -> Result<Means, Error> {
    let mut means = Means::new();
    let mut seen: HashSet<(String, i64, i64, String)> = HashSet::new();
    each_line(path, |line| {
        if strip(line).is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            return Err(Error(format!(
                "{path}: expected chrom, start, end, name, depth sum; got {}",
                py_repr(strip(line))
            )));
        }
        let (chrom, start, end, name) = (c[0], py_int(c[1])?, py_int(c[2])?, c[3]);
        if end <= start {
            return Err(Error(format!(
                "{path}: region {name} {chrom}:{start}-{end} is empty or reversed"
            )));
        }
        if !seen.insert((chrom.to_owned(), start, end, name.to_owned())) {
            return Err(Error(format!(
                "{path}: region {name} {chrom}:{start}-{end} appears twice"
            )));
        }
        means.insert(name.to_owned(), py_float(c[4])? / (end - start) as f64);
        Ok(())
    })?;
    let want: HashSet<(String, i64, i64, String)> = row
        .regions()
        .into_iter()
        .map(|(c, s, e, n)| (c.to_owned(), s, e, n.to_owned()))
        .collect();
    if !seen.is_empty() && seen != want {
        return Err(Error(format!(
            "{path}: its regions differ from the ones `og cn bed --locus {}` prints",
            row.id
        )));
    }
    Ok(means)
}

/// The A, C, G and T reads of one `samtools mpileup` bases column. Read
/// starts (`^` and its MAPQ character), read ends, indels (`+N` and `-N`
/// with their N bases), deletions and skips are not bases; `.` and `,` are
/// the reference base when there is one.
pub fn count_bases(bases: &str, reference: &str) -> [u64; 4] {
    let mut n = [0u64; 4];
    let refb = match reference.chars().next() {
        Some(c) if reference.chars().count() == 1 => c.to_ascii_uppercase(),
        _ => 'N',
    };
    let chars: Vec<char> = bases.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let ch = chars[i];
        if ch == '^' {
            i += 2;
            continue;
        }
        if ch == '+' || ch == '-' {
            let mut j = i + 1;
            while j < chars.len() && chars[j].is_ascii_digit() {
                j += 1;
            }
            let len: usize = chars[i + 1..j]
                .iter()
                .collect::<String>()
                .parse()
                .unwrap_or(0);
            i = j + len;
            continue;
        }
        let b = if ch == '.' || ch == ',' {
            refb
        } else {
            ch.to_ascii_uppercase()
        };
        if let Some(k) = base_index(b) {
            n[k] += 1;
        }
        i += 1;
    }
    n
}

fn base_index(b: char) -> Option<usize> {
    "ACGT".find(b)
}

/// Streams a `samtools mpileup` table (chrom, pos, ref, depth, bases, ...)
/// and hands each line on `chrom` at a position in `want` to `each` with
/// its base counts. A line that is not an mpileup line is an error.
pub fn each_pileup_site(
    path: &str,
    chrom: &str,
    want: &dyn Fn(i64) -> bool,
    mut each: impl FnMut(i64, [u64; 4]),
) -> Result<(), Error> {
    each_line(path, |line| {
        if strip(line).is_empty() {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 || c[1].is_empty() || !c[1].bytes().all(|b| b.is_ascii_digit()) {
            let shown: String = strip(line).chars().take(80).collect();
            return Err(Error(format!(
                "{path}: expected samtools mpileup lines (chrom, pos, ref, depth, bases, ...); got {}",
                py_repr(&shown)
            )));
        }
        let pos = py_int(c[1])?;
        if c[0] == chrom && want(pos) {
            each(pos, count_bases(c[4], c[2]));
        }
        Ok(())
    })
}

/// `{site: (gene-base reads, pseudogene-base reads, other reads)}` of a
/// `call psv` row from `samtools mpileup` over `og cn bed --psv`, each
/// site's two positions pooled: a read that carries the gene's base counts
/// for the gene wherever the aligner put it. A position missing from the
/// pileup has no read.
pub fn read_psv_pileup(path: &str, row: &Locus) -> Result<HashMap<String, SiteCounts>, Error> {
    let mut at: HashMap<i64, usize> = HashMap::new();
    for (i, s) in row.psv.iter().enumerate() {
        at.insert(s.gene_pos + 1, i);
        at.insert(s.pseudo_pos + 1, i);
    }
    let mut counts = vec![[0u64; 4]; row.psv.len()];
    each_pileup_site(path, &row.chrom, &|p| at.contains_key(&p), |pos, n| {
        let k = at[&pos];
        for (a, b) in counts[k].iter_mut().zip(n) {
            *a += b;
        }
    })?;
    Ok(row
        .psv
        .iter()
        .zip(counts)
        .map(|(s, n)| {
            let pick = |bases: &[u8]| -> u64 {
                bases
                    .iter()
                    .filter_map(|b| base_index(*b as char))
                    .map(|k| n[k])
                    .sum()
            };
            let (g, p) = (pick(&s.gene), pick(&s.pseudo));
            (s.name.clone(), (g, p, n.iter().sum::<u64>() - g - p))
        })
        .collect())
}

/// The ALT or HLA contigs named by `samtools view -H` lines, or by contig
/// names one per line: a name ending in `_alt` or starting with `HLA-`
/// (the rule of `scripts/validate-setup.sh`).
pub fn alt_contig(line: &str) -> Option<&str> {
    let name = if let Some(rest) = line.strip_prefix('@') {
        if !rest.starts_with("SQ") {
            return None;
        }
        line.split('\t')
            .find_map(|f| f.strip_prefix("SN:"))
            .unwrap_or("")
    } else {
        strip(line.split('\t').next().unwrap_or(""))
    };
    (name.ends_with("_alt") || name.starts_with("HLA-")).then_some(name)
}

/// The refusal line for a BAM header at `path` (`-` reads stdin) that names
/// an ALT or HLA contig, or `None`.
pub fn check_bam(path: &str) -> Result<Option<String>, Error> {
    let mut found = 0usize;
    let mut first = String::new();
    let mut each = |line: &str| {
        if let Some(name) = alt_contig(line) {
            if found == 0 {
                first = name.to_owned();
            }
            found += 1;
        }
        Ok(())
    };
    if path == "-" {
        crate::each_line_of(std::io::stdin().lock(), "<stdin>", &mut each)?;
    } else {
        each_line(path, &mut each)?;
    }
    Ok((found > 0).then(|| {
        format!(
            "REFUSED: the BAM was aligned to a reference with {found} ALT or HLA contigs (first: {first}). Reads that fit a primary locus and its ALT copy equally get MAPQ 0 or land on the ALT copy, so depth-based copy number reads deletions that are not there. Realign to the default no-ALT reference (docs/realignment.md)."
        )
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bases_skip_starts_ends_and_indels() {
        // ^] read start with its MAPQ; C . , c count; $ ends; +2AG and -1t
        // are indels; * < > N are not bases.
        assert_eq!(count_bases("^]C.,c$+2AGg-1tG*<>N", "C"), [0, 4, 2, 0]);
        assert_eq!(count_bases(".,", "N"), [0, 0, 0, 0]);
    }

    #[test]
    fn alt_and_hla_contigs() {
        assert_eq!(
            alt_contig("@SQ\tSN:chr1_KI270762v1_alt\tLN:1"),
            Some("chr1_KI270762v1_alt")
        );
        assert_eq!(
            alt_contig("HLA-A*01:01:01:01\t3503"),
            Some("HLA-A*01:01:01:01")
        );
        assert_eq!(alt_contig("@SQ\tSN:chr1\tLN:248956422"), None);
        assert_eq!(alt_contig("@HD\tVN:1.6\tSN:x_alt"), None);
    }
}
