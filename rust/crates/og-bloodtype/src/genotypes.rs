//! What the VCF and the reads say at each marker site.
//!
//! The marker table names a few hundred positions; a whole-genome VCF holds
//! millions of records. When the VCF is bgzipped with a `.tbi` or `.csi`
//! beside it (every pipeline VCF is), each position is read through the
//! index: a few hundred BGZF blocks instead of the whole file. Without a
//! usable index (a plain or gzip file, no index, an index older than the
//! VCF, or one that cannot be read) the file is read once, line by line,
//! and a line is split only when its REF covers a marker's position. Both
//! ways apply the same records in the same order, so they give the same
//! sites.
//!
//! Several records can cover one site: an all-sites VCF writes a reference
//! record and then the indel at one position, a split multi-allelic record
//! is two lines, and a deletion that starts upstream spans the site. Every
//! one of them is read, and the answer does not depend on their order:
//!
//! - a record with a non-reference allele that changes the site's base
//!   from upstream makes the site a no call, whatever the site's own
//!   records say;
//! - else, of the records at the site's position, a call of the marker's
//!   ALT wins, then a no call of a record that names the marker's ALT,
//!   then a reference call, then any other no call; the first of equals
//!   counts, and another allele a losing call names stays in the note;
//! - called copies of distinct ALTs listed by one marker add up; repeated
//!   calls of the same ALT count once, and incompatible dosages are no call;
//! - a reference-only record whose REF is shorter than the marker's (an
//!   all-sites record or a gVCF block starting at a deletion marker) says
//!   nothing past its own bases, so the site reads as if it had no record.
//!
//! A call keeps its genotype whatever its quality, but one with a genotype
//! quality under [`MIN_GQ`] or a read depth under [`MIN_DEPTH`] (FORMAT `GQ`,
//! and `DP`, else a gVCF block's `MIN_DP`) is marked weak: the rules then give
//! the systems that use it confidence low. A phased call (`0|1`) keeps the
//! haplotype that carries the marker's ALT, and its phase set (`PS`).

use crate::Error;
use crate::markers::{Marker, Sites};
use indexmap::IndexMap;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// Reads (MAPQ >= 20) before a site missing from the VCF counts as reference,
/// and the read depth (FORMAT `DP`) under which a call is weak.
pub const MIN_DEPTH: f64 = 8.0;

/// The genotype quality (FORMAT `GQ`) under which a call is weak.
pub const MIN_GQ: f64 = 20.0;

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

/// The phase of a call: which of its two haplotypes carries the marker's
/// ALT (`0|1` is `[false, true]`), and its phase set (`PS`), `None` when the
/// record gives none (then the call is phased with the whole chromosome).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Phase {
    pub alt_on: [bool; 2],
    pub block: Option<String>,
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
    /// Why the call is weak (`GQ 5 under 20`), empty when it is not.
    pub weak: String,
    /// The phase of the call, when the VCF gives one.
    pub phase: Option<Phase>,
}

impl Site {
    fn new(row: Marker) -> Site {
        Site {
            row,
            status: Status::Absent,
            n: None,
            depth: None,
            note: String::new(),
            weak: String::new(),
            phase: None,
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
struct ByPos(HashMap<String, BTreeMap<i64, Vec<usize>>>);

impl ByPos {
    fn new(sites: &Genotyped) -> ByPos {
        let mut m: HashMap<String, BTreeMap<i64, Vec<usize>>> = HashMap::new();
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

/// How much a record at a site's own position says about the marker's
/// ALT; the strongest one read counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Weight {
    /// A no call of a record that does not name the marker's ALT.
    Other,
    /// A call with no copy of the marker's ALT.
    Reference,
    /// A no call of a record that names the marker's ALT (`./.`, a
    /// failed FILTER): the ALT may be there.
    Doubt,
    /// A call with the marker's ALT.
    Alt,
}

/// What one record at a site's position says: status, alt copies, note.
#[derive(Debug, Clone, PartialEq)]
struct Reading {
    weight: Weight,
    status: Status,
    n: Option<i64>,
    note: String,
    /// Called marker ALTs and their dosages, independent of VCF allele indices.
    alts: BTreeMap<String, i64>,
    /// Why the call is weak, empty when it is not.
    weak: String,
    phase: Option<Phase>,
}

impl Reading {
    fn no_call(weight: Weight, note: String) -> Reading {
        Reading {
            weight,
            status: Status::NoCall,
            n: None,
            note,
            alts: BTreeMap::new(),
            weak: String::new(),
            phase: None,
        }
    }
}

/// The sample's FORMAT values a call's quality and phase come from.
#[derive(Debug, Clone, Default, PartialEq)]
struct Quality {
    /// `GQ`, as written and as read.
    gq: Option<(String, f64)>,
    /// `DP`, else `MIN_DP`, as written and as read.
    dp: Option<(String, f64)>,
    /// The genotype is written with `|`.
    phased: bool,
    /// `PS`, when given and not `.`.
    ps: Option<String>,
}

impl Quality {
    /// The quality of the sample column `smp` under the keys `fmt`.
    fn of(fmt: &[&str], smp: &[&str], gt: &str) -> Quality {
        let get = |key: &str| {
            fmt.iter()
                .position(|k| *k == key)
                .and_then(|i| smp.get(i))
                .copied()
                .filter(|v| !v.is_empty() && *v != ".")
        };
        let num = |key: &str| {
            get(key).and_then(|v| {
                pycompat::num::parse_float(v)
                    .ok()
                    .filter(|x| x.is_finite())
                    .map(|x| (v.to_owned(), x))
            })
        };
        Quality {
            gq: num("GQ"),
            dp: num("DP").or_else(|| num("MIN_DP")),
            phased: gt.contains('|') && !gt.contains('/'),
            ps: get("PS").map(str::to_owned),
        }
    }

    /// Why a call with this quality is weak, empty when it is not.
    fn weak(&self) -> String {
        let mut why = Vec::new();
        if let Some((raw, gq)) = &self.gq
            && *gq < MIN_GQ
        {
            why.push(format!("GQ {raw} under {}", MIN_GQ as i64));
        }
        if let Some((raw, dp)) = &self.dp
            && *dp < MIN_DEPTH
        {
            why.push(format!("DP {raw} under {}", MIN_DEPTH as i64));
        }
        why.join(", ")
    }
}

/// Everything the records said about one site so far.
#[derive(Debug, Default)]
struct Evidence {
    /// The strongest reading at the site's position.
    best: Option<Reading>,
    /// The notes of losing calls (another allele at the position).
    others: Vec<String>,
    /// The first upstream record that changes the site's base.
    span: Option<String>,
    /// One dosage per called marker ALT, bounded by the marker's ALT list.
    alts: BTreeMap<String, i64>,
    /// Different calls of one ALT disagree about its dosage.
    incompatible: bool,
    /// Why the first weak call of a marker ALT is weak.
    alt_weak: Option<String>,
}

impl Evidence {
    fn add(&mut self, r: Reading) {
        if !r.alts.is_empty() && !r.weak.is_empty() && self.alt_weak.is_none() {
            self.alt_weak = Some(r.weak.clone());
        }
        for (alt, &copies) in &r.alts {
            match self.alts.get(alt) {
                Some(&previous) if previous != copies => self.incompatible = true,
                Some(_) => {}
                None => {
                    self.alts.insert(alt.clone(), copies);
                }
            }
        }
        let loser = match &mut self.best {
            Some(b) if r.weight > b.weight => Some(std::mem::replace(b, r)),
            Some(_) => Some(r),
            None => {
                self.best = Some(r);
                None
            }
        };
        if let Some(l) = loser
            && l.status == Status::Called
            && !l.note.is_empty()
            && !self.others.contains(&l.note)
        {
            self.others.push(l.note);
        }
    }

    /// The site as the records leave it; `Absent` when none said anything.
    fn settle(self, site: &mut Site) {
        if let Some(note) = self.span {
            site.status = Status::NoCall;
            site.n = None;
            site.note = note;
        } else if let Some(b) = self.best {
            let copies: i64 = self.alts.values().sum();
            if self.incompatible || copies > 2 {
                site.status = Status::NoCall;
                site.n = None;
                site.note = "incompatible ALT dosages at this position".to_owned();
                return;
            }
            site.status = b.status;
            site.n = if b.weight == Weight::Alt {
                Some(copies)
            } else {
                b.n
            };
            // The copies of an ALT may come from several records: the call
            // is weak when any of them is, and its phase is the record's only
            // when one ALT was called.
            site.weak = if b.weight == Weight::Alt && b.weak.is_empty() {
                self.alt_weak.clone().unwrap_or_default()
            } else {
                b.weak.clone()
            };
            site.phase = if b.weight == Weight::Alt && self.alts.len() == 1 {
                b.phase.clone()
            } else {
                None
            };
            site.note = if b.status == Status::Called {
                let others = self.others.into_iter().filter(|o| *o != b.note);
                std::iter::once(b.note.clone())
                    .chain(others)
                    .filter(|p| !p.is_empty())
                    .collect::<Vec<_>>()
                    .join("; ")
            } else {
                b.note
            };
        }
    }
}

/// Applies one record line (no line ending) to the sites it covers: those
/// at its position and those its REF spans. `only` keeps the sites at one
/// position (the index path reads each position's records on their own).
fn apply_record(
    line: &str,
    by_pos: &ByPos,
    sites: &Genotyped,
    ev: &mut [Evidence],
    only: Option<i64>,
) {
    let mut c = line.splitn(5, '\t');
    let (Some(chrom), Some(pos), Some(_), Some(rf), Some(_)) =
        (c.next(), c.next(), c.next(), c.next(), c.next())
    else {
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
    let Some(p) = py_pos(pos) else {
        return;
    };
    let span = i64::try_from(rf.len().max(1)).unwrap_or(1);
    let mut hits = here
        .range(p..p.saturating_add(span))
        .filter(|(at, _)| only.is_none_or(|o| o == **at))
        .peekable();
    if hits.peek().is_none() {
        return;
    }
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 10 {
        return;
    }
    let (alts, filt) = (f[4].split(',').collect::<Vec<_>>(), f[6]);
    let fmt: Vec<&str> = f[8].split(':').collect();
    let smp: Vec<&str> = f[9].split(':').collect();
    let gt = fmt
        .iter()
        .position(|k| *k == "GT")
        .and_then(|i| smp.get(i))
        .copied()
        .unwrap_or(".");
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    let quality = Quality::of(&fmt, &smp, gt);
    for (&at, names) in hits {
        if at == p {
            for &i in names {
                if let Some(r) = read_record(&sites[i].row, rf, &alts, filt, &alleles, &quality) {
                    ev[i].add(r);
                }
            }
            continue;
        }
        let offset = usize::try_from(at - p).unwrap_or(0);
        if changes_base(rf, &alts, &alleles, offset) {
            for &i in names {
                ev[i].span.get_or_insert_with(|| {
                    format!(
                        "a record at {}:{} ({rf}>{}, genotype {}) spans this site",
                        f[0],
                        f[1],
                        f[4],
                        alleles.join("/")
                    )
                });
            }
        }
    }
}

/// Whether a record that starts `offset` bases before a site and spans it
/// calls an allele that changes the site's base: one of another length
/// (an indel, a symbolic or `*` allele) or with another base there. A
/// reference or missing genotype changes nothing.
fn changes_base(rf: &str, alts: &[&str], alleles: &[&str], offset: usize) -> bool {
    alleles.iter().any(|a| {
        let Some(k) = a.parse::<usize>().ok().filter(|&k| k > 0) else {
            return false;
        };
        alts.get(k - 1).is_none_or(|alt| {
            alt.len() != rf.len() || alt.as_bytes().get(offset) != rf.as_bytes().get(offset)
        })
    })
}

/// What one record at the site's position says about it; `None` when it
/// says nothing: a reference-only record whose REF is shorter than the
/// marker's, which reads no base past its own.
fn read_record(
    r: &Marker,
    rf: &str,
    alts: &[&str],
    filt: &str,
    alleles: &[&str],
    quality: &Quality,
) -> Option<Reading> {
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
    let doubt = if idx.is_empty() {
        Weight::Other
    } else {
        Weight::Doubt
    };
    let joined = alleles.join("/");
    if alleles.contains(&".") || alleles.len() != 2 {
        return Some(Reading::no_call(doubt, format!("genotype {joined}")));
    }
    let all_ref = alleles.iter().all(|a| *a == "0");
    if rf != r.ref_ {
        if all_ref && rf.len() < r.ref_.len() && r.ref_.starts_with(rf) {
            return None;
        }
        // A record whose REF is not the marker's REF says nothing about
        // the marker's ALT, unless it is a reference call over a span that
        // starts with the marker's REF.
        if !(all_ref && rf.starts_with(r.ref_.as_str())) {
            return Some(Reading::no_call(
                Weight::Other,
                format!("REF {rf} at this position, the marker's REF is {}", r.ref_),
            ));
        }
    }
    if filt != "PASS" && filt != "." && !all_ref {
        return Some(Reading::no_call(doubt, format!("FILTER {filt}")));
    }
    let n = alleles
        .iter()
        .filter(|a| idx.iter().any(|i| i == *a))
        .count() as i64;
    let mut called_alts = BTreeMap::new();
    for (i, alt) in alts.iter().enumerate() {
        if wants.contains(alt) {
            let index = (i + 1).to_string();
            let copies = alleles.iter().filter(|a| **a == index).count() as i64;
            if copies > 0 {
                *called_alts.entry((*alt).to_owned()).or_default() += copies;
            }
        }
    }
    let note = if alleles
        .iter()
        .any(|a| *a != "0" && !idx.iter().any(|i| i == a))
    {
        format!(
            "another allele at this position ({rf}>{}, genotype {joined})",
            alts.join(",")
        )
    } else {
        String::new()
    };
    let phase = quality.phased.then(|| Phase {
        alt_on: [0, 1].map(|h| idx.iter().any(|i| i == alleles[h])),
        block: quality.ps.clone(),
    });
    Some(Reading {
        weight: if n > 0 {
            Weight::Alt
        } else {
            Weight::Reference
        },
        status: Status::Called,
        n: Some(n),
        note,
        alts: called_alts,
        weak: quality.weak(),
        phase,
    })
}

/// Genotypes the sites from the sample's VCF (its one sample column).
pub fn read_vcf(path: &str, markers: Sites) -> Result<(Genotyped, Access), Error> {
    let mut sites: Genotyped = markers
        .into_iter()
        .map(|(k, m)| (k, Site::new(m)))
        .collect();
    let by_pos = ByPos::new(&sites);
    let mut ev: Vec<Evidence> = sites.values().map(|_| Evidence::default()).collect();
    if let Some(records) = indexed_records(path, &sites)? {
        for (pos, line) in &records {
            apply_record(line, &by_pos, &sites, &mut ev, Some(*pos));
        }
        settle(&mut sites, ev);
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
        apply_record(l, &by_pos, &sites, &mut ev, None);
    }
    settle(&mut sites, ev);
    if !header {
        return Err(Error(format!(
            "{path} has no #CHROM header line: not a VCF"
        )));
    }
    Ok((sites, Access::Stream))
}

/// Each site takes what its records said.
fn settle(sites: &mut Genotyped, ev: Vec<Evidence>) {
    for (s, e) in sites.values_mut().zip(ev) {
        e.settle(s);
    }
}

/// The record lines that cover each marker position (start there or span
/// it), read through the VCF's index, with the position they were read
/// for, those of each position in file order; `None` when the file is to
/// be read whole. The `#CHROM` line is checked as a full pass checks it.
fn indexed_records(path: &str, sites: &Genotyped) -> Result<Option<Vec<(i64, String)>>, Error> {
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
                // The index also returns a record whose INFO/END reaches
                // the position; apply_record keeps those whose REF does.
                out.push((*pos, String::from_utf8_lossy(&raw).into_owned()));
            }
        }
    }
    // A site reads only the records that cover its own position, and those
    // come in file order, so the order between positions does not matter.
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

    /// The site as these records, in this order, leave it.
    fn genotype(s: &mut Site, rf: &str, alts: &[&str], filt: &str, alleles: &[&str]) {
        records(s, &[(rf, alts, filt, alleles)]);
    }

    type Rec<'a> = (&'a str, &'a [&'a str], &'a str, &'a [&'a str]);

    fn records(s: &mut Site, recs: &[Rec]) {
        let mut e = Evidence::default();
        for (rf, alts, filt, alleles) in recs {
            if let Some(r) = read_record(&s.row, rf, alts, filt, alleles, &Quality::default()) {
                e.add(r);
            }
        }
        e.settle(s);
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
        genotype(&mut s, "G", &["A"], "PASS", &["0", "1"]);
        assert_eq!(
            (s.status, s.n, s.note.as_str()),
            (Status::Called, Some(1), "")
        );
    }

    #[test]
    fn the_second_alt_of_a_multiallelic_record() {
        let mut s = site("G", "A");
        genotype(&mut s, "G", &["C", "A"], "PASS", &["1", "2"]);
        assert_eq!((s.status, s.n), (Status::Called, Some(1)));
        assert_eq!(
            s.note,
            "another allele at this position (G>C,A, genotype 1/2)"
        );
    }

    #[test]
    fn a_filtered_alt_call_is_no_call_but_a_refcall_is_reference() {
        let mut s = site("G", "A");
        genotype(&mut s, "G", &["A"], "LowQual", &["1", "1"]);
        assert_eq!(
            (s.status, s.note.as_str()),
            (Status::NoCall, "FILTER LowQual")
        );
        let mut s = site("G", "A");
        genotype(&mut s, "G", &["A"], "RefCall", &["0", "0"]);
        assert_eq!((s.status, s.n), (Status::Called, Some(0)));
    }

    #[test]
    fn another_ref_is_no_call_unless_a_reference_span() {
        let mut s = site("G", "A");
        genotype(&mut s, "C", &["T"], "PASS", &["1", "1"]);
        assert_eq!(s.status, Status::NoCall);
        assert_eq!(s.note, "REF C at this position, the marker's REF is G");
        let mut s = site("G", "A");
        genotype(&mut s, "GT", &["G"], "PASS", &["0", "0"]);
        assert_eq!((s.status, s.n), (Status::Called, Some(0)));
    }

    #[test]
    fn a_haploid_or_missing_genotype_is_no_call() {
        for gt in [&["1"][..], &[".", "."], &["0", "1", "1"]] {
            let mut s = site("G", "A");
            genotype(&mut s, "G", &["A"], "PASS", gt);
            assert_eq!(s.status, Status::NoCall, "{gt:?}");
        }
    }

    #[test]
    fn the_strongest_record_at_a_position_counts_in_any_order() {
        let r: Rec = ("T", &["."], ".", &["0", "0"]);
        let ins: Rec = ("T", &["TC"], "PASS", &["0", "1"]);
        let snv: Rec = ("T", &["G"], "PASS", &["1", "0"]);
        let ins_nocall: Rec = ("T", &["TC"], "PASS", &[".", "."]);
        let other_ref: Rec = ("TA", &["T"], "PASS", &["0", "1"]);
        for (recs, want) in [
            (vec![r, ins], (Status::Called, Some(1), "")),
            (vec![ins, r], (Status::Called, Some(1), "")),
            (
                vec![snv, ins],
                (
                    Status::Called,
                    Some(1),
                    "another allele at this position (T>G, genotype 1/0)",
                ),
            ),
            (vec![r, ins_nocall], (Status::NoCall, None, "genotype ./.")),
            (vec![ins_nocall, r], (Status::NoCall, None, "genotype ./.")),
            (vec![other_ref, r], (Status::Called, Some(0), "")),
            (
                vec![other_ref],
                (
                    Status::NoCall,
                    None,
                    "REF TA at this position, the marker's REF is T",
                ),
            ),
            (vec![ins, ins], (Status::Called, Some(1), "")),
        ] {
            let mut s = site("T", "TC");
            records(&mut s, &recs);
            assert_eq!((s.status, s.n, s.note.as_str()), want, "{recs:?}");
        }
    }

    #[test]
    fn a_shorter_reference_record_says_nothing() {
        let mut s = site("CG", "C");
        genotype(&mut s, "C", &["."], ".", &["0", "0"]);
        assert_eq!((s.status, s.n), (Status::Absent, None));
        let mut s = site("CG", "C");
        genotype(&mut s, "C", &["<*>"], ".", &["0", "0"]);
        assert_eq!(s.status, Status::Absent);
        let mut s = site("CG", "C");
        genotype(&mut s, "C", &["T"], "PASS", &["0", "1"]);
        assert_eq!(s.status, Status::NoCall);
    }

    #[test]
    fn a_spanning_allele_changes_the_base_when_its_length_or_base_differs() {
        // CG at the record, the site at offset 1 (the G).
        assert!(changes_base("CG", &["CA"], &["1", "1"], 1));
        assert!(changes_base("CG", &["C"], &["0", "1"], 1));
        assert!(changes_base("CG", &["C", "CA"], &["1", "2"], 1));
        assert!(changes_base("CG", &["*"], &["0", "1"], 1));
        assert!(
            changes_base("CG", &["A"], &["0", "2"], 1),
            "an ALT index the record lacks"
        );
        assert!(!changes_base("CG", &["TG"], &["1", "1"], 1));
        assert!(!changes_base("CG", &["C"], &["0", "0"], 1));
        assert!(!changes_base("CG", &["C"], &[".", "."], 1));
    }

    /// The site as one record with these FORMAT keys and values leaves it.
    fn with_format(r: &str, alt: &str, alleles: &[&str], fmt: &[&str], smp: &[&str]) -> Site {
        let mut s = site(r, alt);
        let mut e = Evidence::default();
        let q = Quality::of(fmt, smp, smp[0]);
        if let Some(rd) = read_record(&s.row, r, &[alt], "PASS", alleles, &q) {
            e.add(rd);
        }
        e.settle(&mut s);
        s
    }

    #[test]
    fn a_call_under_the_gq_or_dp_floor_is_weak_and_keeps_its_genotype() {
        let gt = |g: &str| g.split(['/', '|']).map(str::to_owned).collect::<Vec<_>>();
        for (fmt, smp, want) in [
            (
                &["GT", "GQ", "DP"][..],
                &["0/1", "0", "1"][..],
                "GQ 0 under 20, DP 1 under 8",
            ),
            (&["GT", "GQ", "DP"], &["0/1", "19", "30"], "GQ 19 under 20"),
            (&["GT", "GQ", "DP"], &["0/1", "20", "8"], ""),
            (&["GT", "GQ", "MIN_DP"], &["0/1", "50", "3"], "DP 3 under 8"),
            (&["GT", "DP", "MIN_DP"], &["0/1", "30", "3"], ""),
            (&["GT", "GQ", "DP"], &["0/1", ".", "."], ""),
            (&["GT"], &["0/1"], ""),
        ] {
            let a = gt(smp[0]);
            let alleles: Vec<&str> = a.iter().map(String::as_str).collect();
            let s = with_format("G", "A", &alleles, fmt, smp);
            assert_eq!(
                (s.status, s.n, s.weak.as_str()),
                (Status::Called, Some(1), want),
                "{smp:?}"
            );
        }
        // A weak reference call is weak too; a no call is not marked.
        let s = with_format("G", "A", &["0", "0"], &["GT", "GQ"], &["0/0", "3"]);
        assert_eq!((s.n, s.weak.as_str()), (Some(0), "GQ 3 under 20"));
        let s = with_format("G", "A", &[".", "."], &["GT", "GQ"], &["./.", "3"]);
        assert_eq!((s.status, s.weak.as_str()), (Status::NoCall, ""));
    }

    #[test]
    fn a_phased_call_keeps_the_haplotype_of_the_alt() {
        let s = with_format("G", "A", &["0", "1"], &["GT", "PS"], &["0|1", "100"]);
        assert_eq!(
            s.phase,
            Some(Phase {
                alt_on: [false, true],
                block: Some("100".to_owned())
            })
        );
        let s = with_format("G", "A", &["1", "0"], &["GT", "PS"], &["1|0", "."]);
        assert_eq!(
            s.phase,
            Some(Phase {
                alt_on: [true, false],
                block: None
            })
        );
        let s = with_format("G", "A", &["0", "1"], &["GT"], &["0/1"]);
        assert_eq!(s.phase, None);
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
