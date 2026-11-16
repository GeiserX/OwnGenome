//! The phased reference panel (`panel/chrN.vcf.gz`, one chromosome per
//! file) read into bit-packed haplotypes, for chosen people at chosen
//! sites.
//!
//! Only biallelic SNVs are read: one REF base and one ALT base, each A, C,
//! G or T. Every other record (an indel, a multi-allelic site, no ALT, a
//! `*` allele) is skipped and counted. A chosen person's genotype must be
//! phased, diploid and called (`0|0`, `0|1`, `1|0` or `1|1`), as the
//! panel Beagle phases against is: anything else stops the read, naming
//! the person and the site, so a painting never runs on a guessed allele.
//! People nobody asked for are not read.
//!
//! Haplotype `2 i` is the first allele of person `i` (in the order asked
//! for), haplotype `2 i + 1` the second. A site's alleles are a row of
//! `ceil(2 N / 64)` words, haplotype `h` in bit `h mod 64` of word
//! `h / 64`, 1 for ALT: chromosome 1 of a 470-haplotype panel at 100,000
//! sites is 6 MB.

use crate::{Error, Result, cannot_read, repr};
use memchr::memchr_iter;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// A biallelic SNV.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Site {
    pub pos: u32,
    pub ref_: u8,
    pub alt: u8,
}

impl Site {
    /// The site of a record's POS, REF and ALT when they make a biallelic
    /// SNV.
    pub fn snv(pos: u32, ref_: &[u8], alt: &[u8]) -> Option<Site> {
        let base = |b: &[u8]| match b {
            [c @ (b'A' | b'C' | b'G' | b'T')] => Some(*c),
            _ => None,
        };
        let (r, a) = (base(ref_)?, base(alt)?);
        (r != a).then_some(Site {
            pos,
            ref_: r,
            alt: a,
        })
    }
}

/// The words a row of `haplotypes` bits takes.
pub fn words(haplotypes: usize) -> usize {
    haplotypes.div_ceil(64)
}

/// Bit `h` of a row.
pub fn bit(row: &[u64], h: usize) -> u8 {
    ((row[h / 64] >> (h % 64)) & 1) as u8
}

/// A panel read one site at a time, so a file of any size streams.
pub struct PanelReader {
    path: PathBuf,
    lines: og_vcf::Lines<Box<dyn BufRead + Send>>,
    people: Vec<String>,
    /// (sample column, slot in `people`), by column.
    columns: Vec<(usize, usize)>,
    chrom: Option<String>,
    /// The position read last and the alleles read there.
    at: Option<(u32, Vec<(u8, u8)>)>,
    skipped: u64,
    line: Vec<u8>,
}

impl PanelReader {
    /// Opens `path` (plain, gzip or BGZF) and finds `people` among its
    /// samples. A person missing from the panel, or asked for twice, or
    /// named twice in the panel, is an error.
    pub fn open(path: &Path, people: &[impl AsRef<str>]) -> Result<PanelReader> {
        let input =
            og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(|e| cannot_read(path, &e))?;
        let mut lines = og_vcf::Lines::new(input);
        let header = og_vcf::Header::read(&mut lines)
            .map_err(|e| Error(format!("{}: {e}", path.display())))?;
        let samples: Vec<&[u8]> = header.samples().collect();
        let mut columns = Vec::with_capacity(people.len());
        let mut names = Vec::with_capacity(people.len());
        for (slot, person) in people.iter().enumerate() {
            let person = person.as_ref();
            if names.iter().any(|n| n == person) {
                return Err(Error(format!(
                    "{}: {} is asked for twice",
                    path.display(),
                    repr(person)
                )));
            }
            let mut found = samples
                .iter()
                .enumerate()
                .filter(|(_, s)| **s == person.as_bytes());
            let Some((column, _)) = found.next() else {
                return Err(Error(format!(
                    "{}: no sample {} in the panel",
                    path.display(),
                    repr(person)
                )));
            };
            if found.next().is_some() {
                return Err(Error(format!(
                    "{}: the panel names sample {} twice",
                    path.display(),
                    repr(person)
                )));
            }
            columns.push((column, slot));
            names.push(person.to_owned());
        }
        columns.sort_unstable();
        Ok(PanelReader {
            path: path.to_path_buf(),
            lines,
            people: names,
            columns,
            chrom: None,
            at: None,
            skipped: 0,
            line: Vec::new(),
        })
    }

    /// The people read, in the order asked for.
    pub fn people(&self) -> &[String] {
        &self.people
    }

    /// Two per person.
    pub fn haplotypes(&self) -> usize {
        2 * self.people.len()
    }

    /// The chromosome of the sites read so far.
    pub fn chrom(&self) -> Option<&str> {
        self.chrom.as_deref()
    }

    /// Records skipped so far: not a biallelic SNV, or the same SNV again.
    pub fn skipped(&self) -> u64 {
        self.skipped
    }

    fn fail(&self, what: String) -> Error {
        Error(format!(
            "{}: line {}: {what}",
            self.path.display(),
            self.lines.number()
        ))
    }

    /// The next biallelic SNV, its alleles written into `row` (which holds
    /// [`words`]`(haplotypes)` words); `None` at the end of the file.
    pub fn next_site(&mut self, row: &mut [u64]) -> Result<Option<Site>> {
        assert_eq!(row.len(), words(self.haplotypes()), "row size");
        loop {
            // The line is copied out of the reader, so the reader can say
            // its number while the record is read.
            match self.lines.next_line() {
                Ok(Some(line)) => {
                    self.line.clear();
                    self.line.extend_from_slice(line);
                }
                Ok(None) => return Ok(None),
                Err(e) => return Err(cannot_read(&self.path, &e)),
            }
            if self.line.is_empty() || self.line.starts_with(b"#") {
                continue;
            }
            let rec = og_vcf::Record::parse(&self.line).map_err(|e| {
                Error(format!(
                    "{}: line {}: {e}",
                    self.path.display(),
                    self.lines.number()
                ))
            })?;
            let pos = rec.pos().map_err(|e| {
                Error(format!(
                    "{}: line {}: {e}",
                    self.path.display(),
                    self.lines.number()
                ))
            })?;
            let chrom = String::from_utf8_lossy(rec.chrom());
            match &self.chrom {
                None => self.chrom = Some(chrom.into_owned()),
                Some(c) if *c == chrom => {}
                Some(c) => {
                    let msg = format!(
                        "a record on {chrom} after records on {c}; the panel has one file per chromosome"
                    );
                    return Err(self.fail(msg));
                }
            }
            if let Some((last, _)) = &self.at
                && pos < *last
            {
                let msg = format!("position {pos} after {last}: the panel is not sorted");
                return Err(self.fail(msg));
            }
            let Some(site) = Site::snv(pos, rec.ref_(), rec.alt()) else {
                self.skipped += 1;
                continue;
            };
            match &mut self.at {
                Some((p, seen)) if *p == pos => {
                    if seen.contains(&(site.ref_, site.alt)) {
                        self.skipped += 1;
                        continue;
                    }
                    seen.push((site.ref_, site.alt));
                }
                at => *at = Some((pos, vec![(site.ref_, site.alt)])),
            }
            self.genotypes(&rec, site, row)?;
            return Ok(Some(site));
        }
    }

    /// The chosen people's alleles at `site` into `row`.
    fn genotypes(&self, rec: &og_vcf::Record<'_>, site: Site, row: &mut [u64]) -> Result<()> {
        row.fill(0);
        if self.columns.is_empty() {
            return Ok(());
        }
        if rec.format_keys().next() != Some(b"GT") {
            return Err(self.fail("GT is not the first FORMAT field".into()));
        }
        let Some(samples) = rec.samples() else {
            return Err(self.fail("no sample columns".into()));
        };
        let mut starts = std::iter::once(0).chain(memchr_iter(b'\t', samples).map(|t| t + 1));
        let mut column = 0usize;
        let mut start = starts.next().expect("the first column");
        for &(want, slot) in &self.columns {
            while column < want {
                let Some(next) = starts.next() else {
                    break;
                };
                start = next;
                column += 1;
            }
            if column < want {
                let msg = format!("{} sample columns, fewer than the header names", column + 1);
                return Err(self.fail(msg));
            }
            let field = &samples[start..];
            let end = field
                .iter()
                .position(|&c| c == b':' || c == b'\t')
                .unwrap_or(field.len());
            let gt = &field[..end];
            let (a, b) = match gt {
                [a @ (b'0' | b'1'), b'|', b @ (b'0' | b'1')] => (a - b'0', b - b'0'),
                _ => {
                    let msg = format!(
                        "{} has genotype {} at {}:{}; a panel genotype is phased, diploid and called (0|0, 0|1, 1|0 or 1|1)",
                        self.people[slot],
                        repr(&String::from_utf8_lossy(gt)),
                        self.chrom.as_deref().unwrap_or(""),
                        site.pos
                    );
                    return Err(self.fail(msg));
                }
            };
            for (h, allele) in [(2 * slot, a), (2 * slot + 1, b)] {
                row[h / 64] |= u64::from(allele) << (h % 64);
            }
        }
        Ok(())
    }
}

/// One chromosome of the panel, held whole: the sites and, per site, a row
/// of haplotype bits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Haplotypes {
    chrom: String,
    people: Vec<String>,
    sites: Vec<Site>,
    words: usize,
    bits: Vec<u64>,
}

impl Haplotypes {
    /// Reads the biallelic SNVs of `people` from `path`, every one, or only
    /// those in `keep` (in any order). A site of `keep` the panel lacks is
    /// simply not there: compare [`sites`](Haplotypes::sites) with it.
    pub fn read(
        path: &Path,
        people: &[impl AsRef<str>],
        keep: Option<&[Site]>,
    ) -> Result<Haplotypes> {
        let keep = keep.map(|k| {
            let mut k = k.to_vec();
            k.sort_unstable();
            k
        });
        let mut reader = PanelReader::open(path, people)?;
        let words = words(reader.haplotypes());
        let mut row = vec![0u64; words];
        let (mut sites, mut bits) = (Vec::new(), Vec::new());
        while let Some(site) = reader.next_site(&mut row)? {
            if keep.as_ref().is_none_or(|k| k.binary_search(&site).is_ok()) {
                sites.push(site);
                bits.extend_from_slice(&row);
            }
        }
        Ok(Haplotypes {
            chrom: reader.chrom().unwrap_or("").to_owned(),
            people: reader.people().to_vec(),
            sites,
            words,
            bits,
        })
    }

    /// The chromosome; empty when no record was read.
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn people(&self) -> &[String] {
        &self.people
    }

    pub fn sites(&self) -> &[Site] {
        &self.sites
    }

    /// Two per person.
    pub fn haplotypes(&self) -> usize {
        2 * self.people.len()
    }

    /// The alleles of every haplotype at site `m`.
    pub fn row(&self, m: usize) -> &[u64] {
        &self.bits[m * self.words..(m + 1) * self.words]
    }

    /// The allele (0 REF, 1 ALT) of haplotype `h` at site `m`.
    pub fn allele(&self, m: usize, h: usize) -> u8 {
        assert!(
            h < self.haplotypes(),
            "haplotype {h} of {}",
            self.haplotypes()
        );
        bit(self.row(m), h)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    const HEAD: &str =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    fn panel(d: &Dir, body: &str) -> PathBuf {
        d.put("chr1.vcf.gz", &format!("{HEAD}{body}"))
    }

    #[test]
    fn packs_many_haplotypes_across_words() {
        // 40 people, 80 haplotypes: two words per row.
        let d = Dir::new("panel-wide");
        let names: Vec<String> = (0..40).map(|i| format!("S{i}")).collect();
        let mut text = format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            names.join("\t")
        );
        for pos in 1..=3u32 {
            let gts: Vec<String> = (0..40u32)
                .map(|i| format!("{}|{}", (i + pos) % 2, (i * pos) % 3 % 2))
                .collect();
            text.push_str(&format!(
                "chr1\t{pos}\t.\tA\tC\t.\t.\t.\tGT\t{}\n",
                gts.join("\t")
            ));
        }
        let p = d.put("chr1.vcf", &text);
        let mut asked = names.clone();
        asked.reverse();
        let h = Haplotypes::read(&p, &asked, None).unwrap();
        assert_eq!(h.haplotypes(), 80);
        assert_eq!(h.row(0).len(), 2);
        for m in 0..3 {
            let pos = m as u32 + 1;
            for (slot, name) in asked.iter().enumerate() {
                let i: u32 = name[1..].parse().unwrap();
                assert_eq!(u32::from(h.allele(m, 2 * slot)), (i + pos) % 2);
                assert_eq!(u32::from(h.allele(m, 2 * slot + 1)), (i * pos) % 3 % 2);
            }
        }
    }

    #[test]
    fn keeps_only_the_chosen_sites() {
        let d = Dir::new("panel-keep");
        let p = panel(
            &d,
            "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|0\t1|1\n\
             chr1\t20\t.\tA\tG\t.\t.\t.\tGT\t1|1\t0|0\t1|1\n\
             chr1\t20\t.\tA\tT\t.\t.\t.\tGT\t0|0\t1|0\t0|0\n\
             chr1\t20\t.\tA\tG\t.\t.\t.\tGT\t0|0\t0|0\t0|0\n\
             chr1\t30\t.\tC\tT\t.\t.\t.\tGT\t1|0\t0|1\t0|1\n",
        );
        let all = Haplotypes::read(&p, &["C", "A"], None).unwrap();
        assert_eq!(all.chrom(), "chr1");
        assert_eq!(all.sites().len(), 4, "the second A>G at 20 is skipped");
        let keep = [
            Site::snv(30, b"C", b"T").unwrap(),
            Site::snv(20, b"A", b"T").unwrap(),
            Site::snv(25, b"A", b"T").unwrap(),
        ];
        let h = Haplotypes::read(&p, &["C", "A"], Some(&keep)).unwrap();
        let positions: Vec<(u32, u8)> = h.sites().iter().map(|s| (s.pos, s.alt)).collect();
        assert_eq!(positions, [(20, b'T'), (30, b'T')]);
        // C then A: C is 0|0 at 20 A>T, A is 0|0; at 30 C is 0|1, A 1|0.
        assert_eq!(h.row(0), [0b0000]);
        assert_eq!(h.row(1), [0b0110]);
    }

    #[test]
    fn refuses_what_it_cannot_read_truly() {
        let d = Dir::new("panel-bad");
        let err = |body: &str, people: &[&str]| {
            let p = panel(&d, body);
            Haplotypes::read(&p, people, None).unwrap_err().0
        };
        let ok = "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|0\t1|1\n";
        let e = err(ok, &["D"]);
        assert!(e.ends_with("no sample 'D' in the panel"), "{e}");
        let e = err(ok, &["A", "A"]);
        assert!(e.ends_with("'A' is asked for twice"), "{e}");
        let e = err("chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0|0\t1|1\n", &["A"]);
        assert!(
            e.ends_with("line 3: A has genotype '0/1' at chr1:10; a panel genotype is phased, diploid and called (0|0, 0|1, 1|0 or 1|1)"),
            "{e}"
        );
        let e = err("chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t.|.\t1|1\n", &["B"]);
        assert!(e.contains("B has genotype '.|.'"), "{e}");
        let e = err("chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1\t1|1\n", &["B"]);
        assert!(e.contains("B has genotype '1'"), "{e}");
        let e = err("chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|0\n", &["C"]);
        assert!(
            e.contains("2 sample columns, fewer than the header names"),
            "{e}"
        );
        let e = err(
            "chr1\t10\t.\tA\tG\t.\t.\t.\tDP:GT\t3:0|1\t3:0|0\t3:1|1\n",
            &["C"],
        );
        assert!(e.contains("GT is not the first FORMAT field"), "{e}");
        let e = err(
            &format!("{ok}chr2\t5\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|0\t1|1\n"),
            &["A"],
        );
        assert!(
            e.contains("line 4: a record on chr2 after records on chr1"),
            "{e}"
        );
        let e = err(
            &format!("{ok}chr1\t9\t.\tA\tG\t.\t.\t.\tGT\t0|1\t0|0\t1|1\n"),
            &["A"],
        );
        assert!(
            e.contains("position 9 after 10: the panel is not sorted"),
            "{e}"
        );
        // A person not asked for may be anything.
        let p = panel(&d, "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t./.\t1|1\n");
        assert_eq!(
            Haplotypes::read(&p, &["A", "C"], None).unwrap().row(0),
            [0b1110]
        );
        // No one asked for: the sites alone.
        let h = Haplotypes::read(&p, &[] as &[&str], None).unwrap();
        assert_eq!((h.sites().len(), h.haplotypes()), (1, 0));
    }

    #[test]
    fn snv_rule() {
        assert!(Site::snv(1, b"A", b"G").is_some());
        for (r, a) in [
            (&b"A"[..], &b"A"[..]),
            (b"AT", b"A"),
            (b"A", b"G,T"),
            (b"A", b"."),
            (b"A", b"*"),
            (b"a", b"g"),
            (b"N", b"A"),
        ] {
            assert!(Site::snv(1, r, a).is_none(), "{r:?} {a:?}");
        }
    }
}
