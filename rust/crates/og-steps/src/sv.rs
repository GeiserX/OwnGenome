//! `og sv add-callers` and `og sv postfilter`: step 05's structural
//! variants on their way into AnnotSV and the checks a class 4 or 5 call
//! needs after it (port of `bin/cnvpytor_to_vcf.py` and
//! `bin/annotsv_postfilter.py`).
//!
//! - [`add_callers`]: the base VCF (Manta's calls) again, plus the Delly
//!   and CNVpytor deletions and duplications no call already in it covers.
//! - [`postfilter`]: AnnotSV's class 4 and 5 rows with their caller, the
//!   exon-edge, population and depth checks, and the region rules.
//!
//! This module holds the rules both share, and that step 39's couple
//! screen (`bin/couple_report.py`) applies too: the read-depth rule for a
//! long Delly call ([`depth_contradicts`]), whether a record is a call the
//! sample carries ([`call_holds`]) and the reciprocal-overlap index
//! ([`OverlapIndex`]).

pub mod add_callers;
pub mod postfilter;

use std::collections::{BTreeMap, HashMap};

/// A Delly DEL or DUP longer than this many bp needs its read depth to
/// agree: its FORMAT RDCN (Delly's read-depth copy number) must be below
/// (DEL) or above (DUP) the copies a normal region has. 1,000 bp is
/// Delly's own line between small and large SVs in its read counts
/// (`src/delly.h`, `indelsize`); short-read Delly can mark PASS a long SV
/// whose RDCN shows no change, and a short one rests on split reads alone.
pub const DEPTH_MIN_SIZE: i64 = 1000;

/// The type of a call, as the rules read it: `DUP:TANDEM` is a DUP.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SvKind {
    Del,
    Dup,
    Other,
}

impl SvKind {
    /// The symbolic ALT the type writes (`<DEL>`, `<DUP>`).
    pub fn alt(self) -> &'static str {
        match self {
            SvKind::Del => "<DEL>",
            SvKind::Dup => "<DUP>",
            SvKind::Other => ".",
        }
    }
}

/// The sample's sex, which sets the normal copies on chrX and chrY.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    /// Not given: 1 and 2 copies both count as normal on chrX and chrY
    /// outside the PARs.
    Unknown,
}

impl Sex {
    /// `male`, `female` or `` (not given); anything else is `None`.
    pub fn parse(s: &str) -> Option<Sex> {
        match s {
            "male" => Some(Sex::Male),
            "female" => Some(Sex::Female),
            "" => Some(Sex::Unknown),
            _ => None,
        }
    }
}

/// Delly's read-depth copy number of a call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rdcn {
    /// The record has no RDCN field (not Delly 2.x): never held against
    /// the call.
    Absent,
    /// The field is there but says nothing: `.`, empty, negative (no reads
    /// in the flanks) or not a whole number of copies.
    NoEvidence,
    /// A whole number of copies.
    Copies(u64),
}

impl Rdcn {
    /// RDCN of a Delly record's first sample: `keys` are FORMAT's, `vals`
    /// the sample's. The first `RDCN` key counts, and a value the sample
    /// column does not reach is `.`.
    pub fn from_format(keys: &[&str], vals: &[&str]) -> Rdcn {
        let Some(i) = keys.iter().position(|&k| k == "RDCN") else {
            return Rdcn::Absent;
        };
        let cn = vals.get(i).copied().unwrap_or(".");
        if cn.is_empty() || !cn.bytes().all(|b| b.is_ascii_digit()) {
            return Rdcn::NoEvidence;
        }
        // Digits only: a value past u64 is more copies than any normal.
        Rdcn::Copies(cn.parse().unwrap_or(u64::MAX))
    }

    /// RDCN read as a number (an AnnotSV cell): a finite whole number of
    /// zero or more is copies, any other number no evidence.
    pub fn from_number(v: f64) -> Rdcn {
        if v.is_finite() && v >= 0.0 && v.fract() == 0.0 {
            // `as` saturates: 1e300 copies is u64::MAX copies.
            Rdcn::Copies(v as u64)
        } else {
            Rdcn::NoEvidence
        }
    }
}

/// Whether a Delly DEL or DUP of `len` bp has a read depth that does not
/// show the change: not below the normal copies for a DEL, not above them
/// for a DUP. The normal is 1 on a male's chrX or chrY outside the PARs and
/// 2 elsewhere; with no sex, either 1 or 2 there, so only an RDCN at or
/// past both is held against the call. A call of [`DEPTH_MIN_SIZE`] bp or
/// less, or with no RDCN field, is never contradicted; an RDCN that says
/// nothing ([`Rdcn::NoEvidence`]) always is.
pub fn depth_contradicts(
    kind: SvKind,
    chrom: &str,
    len: i64,
    rdcn: Rdcn,
    sex: Sex,
    in_par: bool,
) -> bool {
    if len <= DEPTH_MIN_SIZE {
        return false;
    }
    let cn = match rdcn {
        Rdcn::Absent => return false,
        Rdcn::NoEvidence => return true,
        Rdcn::Copies(n) => n,
    };
    let sex_chrom = matches!(chrom, "chrX" | "chrY") && !in_par;
    let (low, high) = match (sex_chrom, sex) {
        (true, Sex::Male) => (1, 1),
        (true, Sex::Unknown) => (1, 2),
        _ => (2, 2),
    };
    if kind == SvKind::Del {
        cn >= high
    } else {
        cn <= low
    }
}

/// The GT of a record's first sample (`./.` when there is none).
pub fn first_gt<'a>(fields: &[&'a str]) -> &'a str {
    if fields.len() < 10 {
        return "./.";
    }
    let Some(i) = fields[8].split(':').position(|k| k == "GT") else {
        return "./.";
    };
    match fields[9].split(':').nth(i) {
        Some(v) if !v.is_empty() => v,
        _ => "./.",
    }
}

/// Whether a record is a call this sample carries: FILTER `PASS` or `.`,
/// and in the first sample an allele other than `0` or `.`. A record
/// Manta filtered or genotyped as reference or no-call is not; a
/// sites-only record (`n_samples` 0) has no genotype and is.
pub fn call_holds(fields: &[&str], n_samples: usize) -> bool {
    if !matches!(fields.get(6).copied(), Some("PASS" | ".")) {
        return false;
    }
    if n_samples == 0 {
        return true;
    }
    first_gt(fields)
        .split(['/', '|'])
        .any(|a| !matches!(a, "0" | "." | ""))
}

/// Calls by contig and type, for the reciprocal-overlap check: each group
/// keeps its calls sorted by start, with the longest call's length, so a
/// query looks only at the calls that can reach it.
#[derive(Debug, Default)]
pub struct OverlapIndex {
    groups: HashMap<(String, SvKind), Group>,
}

#[derive(Debug, Default)]
struct Group {
    by_start: BTreeMap<i64, Vec<i64>>,
    max_len: i64,
}

impl OverlapIndex {
    /// Adds the call `start..end`.
    pub fn add(&mut self, chrom: &str, kind: SvKind, start: i64, end: i64) {
        let g = self.groups.entry((chrom.to_owned(), kind)).or_default();
        g.by_start.entry(start).or_default().push(end);
        g.max_len = g.max_len.max(end - start);
    }

    /// Whether a call of the same type on the same contig overlaps
    /// `start..end` by at least `frac` of the longer of the two.
    pub fn covers(&self, chrom: &str, kind: SvKind, start: i64, end: i64, frac: f64) -> bool {
        let Some(g) = self.groups.get(&(chrom.to_owned(), kind)) else {
            return false;
        };
        let len = (end - start).max(1);
        let lo = start - g.max_len;
        if lo >= end {
            return false;
        }
        g.by_start.range(lo..end).any(|(&s, ends)| {
            ends.iter().any(|&e| {
                let overlap = e.min(end) - s.max(start);
                overlap > 0 && overlap as f64 >= frac * (e - s).max(len) as f64
            })
        })
    }
}

/// A chromosome name without `chr` (any case): AnnotSV writes them so,
/// step 22 with it.
pub fn chrom_key(c: &str) -> &str {
    match c.get(..3) {
        Some(p) if p.eq_ignore_ascii_case("chr") => &c[3..],
        _ => c,
    }
}

/// The value of the first `key=value` item of a `;`-separated INFO, or
/// `""` (also for a flag).
pub fn info_value<'a>(info: &'a str, key: &str) -> &'a str {
    info.split(';')
        .find_map(|kv| {
            let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
            (k == key).then_some(v)
        })
        .unwrap_or("")
}

#[cfg(test)]
mod tests;
