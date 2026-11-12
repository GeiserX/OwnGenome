//! What the VCF and the reads say at each marker site.
//!
//! The marker table names a few hundred positions; a whole-genome VCF holds
//! millions of records. When the VCF is bgzipped with a `.tbi` or `.csi`
//! beside it (every pipeline VCF is), each position is read through the
//! index: a few hundred BGZF blocks instead of the whole file. Without a
//! usable index (a plain or gzip file, no index, an index older than the
//! VCF, or one that cannot be read) the file is read once, line by line,
//! and a line is split only when its CHROM and POS are a marker's. Both
//! ways apply the same records in the same order, so they give the same
//! sites.

use crate::Error;
use crate::markers::{Marker, Sites};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::path::Path;

/// Reads (MAPQ >= 20) before a site missing from the VCF counts as reference.
pub const MIN_DEPTH: f64 = 8.0;

/// What is known at one site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// No VCF record yet.
    Absent,
    /// A record with a genotype the site can read.
    Called,
    /// A record that says nothing about the site's alleles.
    NoCall,
    /// No record and no BAM: read as reference.
    Assumed,
    /// No record, and enough reads: reference.
    Covered,
    /// No record and too few reads.
    LowDepth,
}

impl Status {
    /// The word the evidence table and the notes use.
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Absent => "absent",
            Status::Called => "called",
            Status::NoCall => "no call",
            Status::Assumed => "assumed",
            Status::Covered => "covered",
            Status::LowDepth => "low depth",
        }
    }
}

/// What the VCF and the reads say at one site: status, alt copies (`None`
/// when not known), read depth, note.
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub row: Marker,
    pub status: Status,
    pub n: Option<i64>,
    pub depth: Option<f64>,
    pub note: String,
}

impl Site {
    fn new(row: Marker) -> Site {
        Site {
            row,
            status: Status::Absent,
            n: None,
            depth: None,
            note: String::new(),
        }
    }
}

/// The sites by name, in table order.
pub type Genotyped = IndexMap<String, Site>;

/// How the VCF was read, for the log and the tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Access {
    /// Through its `.tbi` or `.csi`.
    Index,
    /// One pass over every line.
    Stream,
}

/// `chr1` stays, `1` becomes `chr1`, `MT` and `M` become `chrM`.
pub fn norm_chrom(c: &str) -> String {
    if c.starts_with("chr") {
        c.to_owned()
    } else if c == "MT" || c == "M" {
        "chrM".to_owned()
    } else {
        format!("chr{c}")
    }
}

/// The sites at each (chrom, pos), in table order.
struct ByPos(HashMap<String, HashMap<i64, Vec<usize>>>);

impl ByPos {
    fn new(sites: &Genotyped) -> ByPos {
        let mut m: HashMap<String, HashMap<i64, Vec<usize>>> = HashMap::new();
        for (i, s) in sites.values().enumerate() {
            m.entry(s.row.chrom.clone())
                .or_default()
                .entry(s.row.pos)
                .or_default()
                .push(i);
        }
        ByPos(m)
    }
}

/// `int(s)` as Python reads a POS.
fn py_pos(s: &str) -> Option<i64> {
    pycompat::num::parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
}

/// The sample count of a `#CHROM` line: an error unless it is one.
fn check_header(line: &str, path: &str) -> Result<(), Error> {
    let n = line.split('\t').count() as i64 - 9;
    if n < 1 {
        return Err(Error(format!(
            "{path} has no sample column: a genotype is needed at every site"
        )));
    }
    if n > 1 {
        return Err(Error(format!(
            "{path} has {n} samples: this step reads one person's VCF"
        )));
    }
    Ok(())
}

/// Applies one record line (no line ending) to the sites at its position.
fn apply_record(line: &str, by_pos: &ByPos, sites: &mut Genotyped) {
    let mut c = line.splitn(3, '\t');
    let (Some(chrom), Some(pos), Some(_)) = (c.next(), c.next(), c.next()) else {
        return;
    };
    let chrom = if chrom.starts_with("chr") {
        std::borrow::Cow::Borrowed(chrom)
    } else {
        std::borrow::Cow::Owned(norm_chrom(chrom))
    };
    let Some(here) = by_pos.0.get(chrom.as_ref()) else {
        return;
    };
    let Some(names) = py_pos(pos).and_then(|p| here.get(&p)) else {
        return;
    };
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 10 {
        return;
    }
    let (rf, alts, filt) = (f[3], f[4].split(',').collect::<Vec<_>>(), f[6]);
    let fmt: Vec<&str> = f[8].split(':').collect();
    let smp: Vec<&str> = f[9].split(':').collect();
    let gt = fmt
        .iter()
        .position(|k| *k == "GT")
        .and_then(|i| smp.get(i))
        .copied()
        .unwrap_or(".");
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    for &i in names {
        let site = &mut sites[i];
        if site.status == Status::Called {
            continue; // a second record at the position: the first one counts
        }
        genotype_site(site, rf, &alts, filt, &alleles);
    }
}

/// The site's status, alt copies and note from one record.
pub fn genotype_site(site: &mut Site, rf: &str, alts: &[&str], filt: &str, alleles: &[&str]) {
    let r = &site.row;
    // A null row may list several ALTs (null alleles at one position): any of them counts.
    let wants: Vec<&str> = if r.ref_ == rf {
        r.alt.split(',').collect()
    } else {
        Vec::new()
    };
    let idx: Vec<String> = alts
        .iter()
        .enumerate()
        .filter(|(_, a)| wants.contains(a))
        .map(|(i, _)| (i + 1).to_string())
        .collect();
    let joined = alleles.join("/");
    if alleles.contains(&".") || alleles.len() != 2 {
        site.status = Status::NoCall;
        site.note = format!("genotype {joined}");
        return;
    }
    // A record whose REF is not the marker's REF says nothing about the
    // marker's ALT, unless it is a reference call over a span that starts
    // with the marker's REF.
    if rf != r.ref_ && !(rf.starts_with(r.ref_.as_str()) && alleles.iter().all(|a| *a == "0")) {
        site.note = format!("REF {rf} at this position, the marker's REF is {}", r.ref_);
        site.status = Status::NoCall;
        return;
    }
    if filt != "PASS" && filt != "." && alleles.iter().any(|a| *a != "0") {
        site.status = Status::NoCall;
        site.note = format!("FILTER {filt}");
        return;
    }
    site.status = Status::Called;
    site.n = Some(
        alleles
            .iter()
            .filter(|a| idx.iter().any(|i| i == *a))
            .count() as i64,
    );
    site.note.clear();
    if alleles
        .iter()
        .any(|a| *a != "0" && !idx.iter().any(|i| i == a))
    {
        site.note = format!(
            "another allele at this position ({rf}>{}, genotype {joined})",
            alts.join(",")
        );
    }
}

/// Genotypes the sites from the sample's VCF (its one sample column).
pub fn read_vcf(path: &str, markers: Sites) -> Result<(Genotyped, Access), Error> {
    let mut sites: Genotyped = markers
        .into_iter()
        .map(|(k, m)| (k, Site::new(m)))
        .collect();
    let by_pos = ByPos::new(&sites);
    if let Some(records) = indexed_records(path, &sites)? {
        for line in &records {
            apply_record(line, &by_pos, &mut sites);
        }
        return Ok((sites, Access::Index));
    }
    let mut lines = og_io::Lines::new(
        og_io::open(path).map_err(|e| Error::io(&e, path))?,
        og_io::Utf8::Replace,
    );
    let mut header = false;
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        let l = line.strip_suffix('\n').unwrap_or(&line);
        if l.starts_with("##") {
            continue;
        }
        if l.starts_with("#CHROM") {
            check_header(l, path)?;
            header = true;
            continue;
        }
        apply_record(l, &by_pos, &mut sites);
    }
    if !header {
        return Err(Error(format!(
            "{path} has no #CHROM header line: not a VCF"
        )));
    }
    Ok((sites, Access::Stream))
}

/// The record lines at the marker positions, read through the VCF's
/// index, those of each position in file order; `None` when the file is to
/// be read whole. The `#CHROM` line is checked as a full pass checks it.
fn indexed_records(path: &str, sites: &Genotyped) -> Result<Option<Vec<String>>, Error> {
    let Some(mut reader) = usable_index(path) else {
        return Ok(None);
    };
    // The header comes first in a bgzipped VCF: read it up to #CHROM. A
    // record before any #CHROM line is a file a full pass reads its own way.
    let mut lines = og_io::Lines::new(
        og_io::open(path).map_err(|e| Error::io(&e, path))?,
        og_io::Utf8::Replace,
    );
    let mut line = String::new();
    loop {
        if !lines
            .read_line(&mut line)
            .map_err(|e| Error::io(&e, path))?
        {
            return Ok(None);
        }
        let l = line.strip_suffix('\n').unwrap_or(&line);
        if l.starts_with("#CHROM") {
            check_header(l, path)?;
            break;
        }
        if !l.starts_with("##") {
            return Ok(None);
        }
    }
    drop(lines);
    // Each marker position once, in table order; the index's names that
    // read as its chromosome, in file order.
    let mut wanted: IndexMap<(String, i64), ()> = IndexMap::new();
    for s in sites.values() {
        wanted.insert((s.row.chrom.clone(), s.row.pos), ());
    }
    let names: Vec<String> = reader.index().names().to_vec();
    let mut out = Vec::new();
    for (chrom, pos) in wanted.keys() {
        let Ok(beg) = u64::try_from(*pos - 1) else {
            continue;
        };
        for name in names.iter().filter(|n| norm_chrom(n) == *chrom) {
            for got in reader.query(name, beg, beg + 1) {
                let Ok(raw) = got else {
                    return Ok(None);
                };
                let l = String::from_utf8_lossy(&raw).into_owned();
                if l.split('\t').nth(1).and_then(py_pos) == Some(*pos) {
                    out.push(l);
                }
            }
        }
    }
    // A site reads only the records at its own position, and those come in
    // file order, so the order between positions does not matter.
    Ok(Some(out))
}

/// The VCF's index reader, when there is one to trust: a `.tbi` or `.csi`
/// beside it, for VCF lines, no older than the VCF.
fn usable_index(path: &str) -> Option<og_io::tabix::Reader> {
    let vcf_time = std::fs::metadata(path).ok()?.modified().ok()?;
    let index = [".tbi", ".csi"].into_iter().find_map(|ext| {
        let p = format!("{path}{ext}");
        Path::new(&p).is_file().then_some(p)
    })?;
    let index_time = std::fs::metadata(&index).ok()?.modified().ok()?;
    if index_time < vcf_time {
        return None;
    }
    let reader = og_io::tabix::Reader::with_index(path, &index).ok()?;
    (reader.index().conf().preset == og_io::tabix::Preset::Vcf).then_some(reader)
}

/// `samtools bedcov` output over the marker sites (name in column 6) and the
/// RHD regions (a 4-column BED, name in column 4): mean depth by name. A
/// line with fewer than 5 columns, a `#` line, or one whose numbers do not
/// read is skipped.
pub fn read_depth(path: &str) -> Result<HashMap<String, f64>, Error> {
    let mut depth = HashMap::new();
    let file = std::fs::File::open(path).map_err(|e| Error::io(&e, path))?;
    let lines = og_io::Lines::new(std::io::BufReader::new(file), og_io::Utf8::Strict);
    for l in lines {
        let l = l.map_err(|e| Error::io(&e, path))?;
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 5 || c[0].starts_with('#') {
            continue;
        }
        let (Some(start), Some(end), Ok(total)) = (
            py_pos(c[1]),
            py_pos(c[2]),
            pycompat::num::parse_float(c[c.len() - 1]),
        ) else {
            continue;
        };
        let name = if c.len() >= 7 { c[5] } else { c[3] };
        depth.insert(name.to_owned(), total / (end - start).max(1) as f64);
    }
    Ok(depth)
}

/// A site with no VCF record: reference when the reads cover it, else not
/// known. Every site gets its depth.
pub fn apply_depth(sites: &mut Genotyped, depth: Option<&HashMap<String, f64>>) {
    for s in sites.values_mut() {
        let d = depth.and_then(|m| m.get(&s.row.name).copied());
        s.depth = d;
        if s.status != Status::Absent {
            continue;
        }
        match (depth, d) {
            (None, _) => {
                s.status = Status::Assumed;
                s.n = Some(0);
                s.note = "no VCF record; no BAM depth to confirm the reference".to_owned();
            }
            (Some(_), None) => {
                s.status = Status::NoCall;
                s.note = "no VCF record and no depth row (the BAM has no such contig)".to_owned();
            }
            (Some(_), Some(d)) if d >= MIN_DEPTH => {
                s.status = Status::Covered;
                s.n = Some(0);
            }
            (Some(_), Some(d)) => {
                s.status = Status::LowDepth;
                s.note = format!(
                    "no VCF record and {} reads at MAPQ 20 (fewer than {})",
                    pycompat::fmt::fixed(d, 0),
                    MIN_DEPTH as i64
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(r: &str, alt: &str) -> Site {
        Site::new(Marker {
            ref_: r.into(),
            alt: alt.into(),
            ..Marker::default()
        })
    }

    #[test]
    fn chroms_normalise() {
        assert_eq!(norm_chrom("chr7"), "chr7");
        assert_eq!(norm_chrom("7"), "chr7");
        assert_eq!(norm_chrom("MT"), "chrM");
        assert_eq!(norm_chrom("M"), "chrM");
    }

    #[test]
    fn a_het_call_counts_one_alt_copy() {
        let mut s = site("G", "A");
        genotype_site(&mut s, "G", &["A"], "PASS", &["0", "1"]);
        assert_eq!(
            (s.status, s.n, s.note.as_str()),
            (Status::Called, Some(1), "")
        );
    }

    #[test]
    fn the_second_alt_of_a_multiallelic_record() {
        let mut s = site("G", "A");
        genotype_site(&mut s, "G", &["C", "A"], "PASS", &["1", "2"]);
        assert_eq!((s.status, s.n), (Status::Called, Some(1)));
        assert_eq!(
            s.note,
            "another allele at this position (G>C,A, genotype 1/2)"
        );
    }

    #[test]
    fn a_filtered_alt_call_is_no_call_but_a_refcall_is_reference() {
        let mut s = site("G", "A");
        genotype_site(&mut s, "G", &["A"], "LowQual", &["1", "1"]);
        assert_eq!(
            (s.status, s.note.as_str()),
            (Status::NoCall, "FILTER LowQual")
        );
        let mut s = site("G", "A");
        genotype_site(&mut s, "G", &["A"], "RefCall", &["0", "0"]);
        assert_eq!((s.status, s.n), (Status::Called, Some(0)));
    }

    #[test]
    fn another_ref_is_no_call_unless_a_reference_span() {
        let mut s = site("G", "A");
        genotype_site(&mut s, "C", &["T"], "PASS", &["1", "1"]);
        assert_eq!(s.status, Status::NoCall);
        assert_eq!(s.note, "REF C at this position, the marker's REF is G");
        let mut s = site("G", "A");
        genotype_site(&mut s, "GT", &["G"], "PASS", &["0", "0"]);
        assert_eq!((s.status, s.n), (Status::Called, Some(0)));
    }

    #[test]
    fn a_haploid_or_missing_genotype_is_no_call() {
        for gt in [&["1"][..], &[".", "."], &["0", "1", "1"]] {
            let mut s = site("G", "A");
            genotype_site(&mut s, "G", &["A"], "PASS", gt);
            assert_eq!(s.status, Status::NoCall, "{gt:?}");
        }
    }

    #[test]
    fn depth_decides_a_site_with_no_record() {
        let mut g: Genotyped = IndexMap::new();
        for (name, _) in [("a", 0), ("b", 0), ("c", 0)] {
            let mut s = site("G", "A");
            s.row.name = name.into();
            g.insert(name.into(), s);
        }
        let depth: HashMap<String, f64> = [("a".to_owned(), 8.0), ("b".to_owned(), 7.9)].into();
        apply_depth(&mut g, Some(&depth));
        assert_eq!(g["a"].status, Status::Covered);
        assert_eq!(g["b"].status, Status::LowDepth);
        assert_eq!(
            g["b"].note,
            "no VCF record and 8 reads at MAPQ 20 (fewer than 8)"
        );
        assert_eq!(g["c"].status, Status::NoCall);
        apply_depth(&mut g, None);
        assert_eq!((g["a"].status, g["a"].depth), (Status::Covered, None));
    }
}
