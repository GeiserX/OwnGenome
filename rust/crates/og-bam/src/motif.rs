//! `og str motif`: which repeat motif do the reads at RFC1 (or another
//! motif locus) carry? The port of `bin/rfc1_motif.py`.
//!
//! ExpansionHunter sizes RFC1 with the degenerate unit AARRG, so a common
//! benign expansion (AAAAG, AAAGG) is flagged as readily as a pathogenic one
//! (AAGGG, CANVAS, Cortese et al. 2019; ACAGG). Short reads cannot size a
//! long expansion, but they can say which motif it is made of.
//!
//! For every read at the locus this takes the read's bases over the repeat
//! and [`FLANK`] bases either side, with the soft clip next to them, and the
//! whole sequence of an unmapped read placed within [`FETCH_PAD`] bases and
//! of a mate mapped elsewhere (mate rescue). The longest tandem run of one
//! unit in those bases counts for its motif when it is longer than the
//! reference allele. A motif is named by its known form whatever its
//! rotation and strand. Duplicate, secondary, supplementary and QC-failed
//! alignments do not count.
//!
//! ```text
//! og str motif --vcf STRANGER_VCF --out OUT.tsv (--bam BAM | --skip REASON)
//!              [--catalog CATALOG.json [--flagged]] [--no-softclip] [--no-mates]
//! ```
//!
//! The output is the Python's, byte for byte: a header and one row per
//! locus (`locus chrom start end repeat_reads dominant_motif dominant_class
//! dominant_share motif_shares verdict note`).

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

use indexmap::IndexMap;
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::fmt::{fixed, percent};
use pycompat::json::{self, Value};
use pycompat::num::parse_int;
use pycompat::text::{is_decimal, isdigit};
use pycompat::value::{py_str, truthy};

use crate::Bam;
use crate::args::Args;
use rust_htslib::bam::Record;

/// Reference bases taken either side of the repeat.
pub const FLANK: i64 = 50;
/// The fetched window around the repeat: unmapped reads placed there count.
pub const FETCH_PAD: i64 = 1000;
const MIN_UNITS: i64 = 4;
const MIN_READS: usize = 3;
const PATHOGENIC_SHARE: f64 = 0.2;
const BENIGN_SHARE: f64 = 0.5;

const COLUMNS: &str = "locus\tchrom\tstart\tend\trepeat_reads\tdominant_motif\tdominant_class\t\
                       dominant_share\tmotif_shares\tverdict\tnote";

/// What a known motif means at its locus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    Benign,
    Pathogenic,
}

impl Class {
    fn name(self) -> &'static str {
        match self {
            Class::Benign => "benign",
            Class::Pathogenic => "pathogenic",
        }
    }
}

/// How one locus is read.
#[derive(Debug, Clone)]
pub struct Spec {
    pub locus: String,
    /// Known motifs in catalog order (benign first, a motif in both lists
    /// keeps its first place and reads pathogenic, as a Python dict does).
    pub known: IndexMap<String, Class>,
    pub period: usize,
    pub disease: String,
    pub recessive: bool,
    /// The catalog's span (contig, 0-based start, end); `None` takes the VCF's.
    pub span: Option<(String, i64, i64)>,
}

impl Spec {
    /// RFC1 read without a catalog.
    pub fn rfc1() -> Self {
        Spec {
            locus: "RFC1".into(),
            known: [
                ("AAAAG", Class::Benign),
                ("AAAGG", Class::Benign),
                ("AAGGG", Class::Pathogenic),
                ("ACAGG", Class::Pathogenic),
            ]
            .into_iter()
            .map(|(m, c)| (m.to_owned(), c))
            .collect(),
            period: 5,
            disease: "CANVAS".into(),
            recessive: true,
            span: None,
        }
    }

    /// Every form of each known motif mapped to its name; a later motif
    /// wins a form two share.
    fn names(&self) -> HashMap<Vec<u8>, String> {
        let mut out = HashMap::new();
        for k in self.known.keys() {
            for f in forms(k.as_bytes()) {
                out.insert(f, k.clone());
            }
        }
        out
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Every rotation of a motif on both strands, without repeats.
pub fn forms(m: &[u8]) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    for x in [m.to_vec(), revcomp(m)] {
        for i in 0..x.len().max(1) {
            let r = [&x[i.min(x.len())..], &x[..i.min(x.len())]].concat();
            if !out.contains(&r) {
                out.push(r);
            }
        }
    }
    out
}

/// The known name of a motif, or its smallest form on either strand.
fn motif_name(m: &[u8], names: &HashMap<Vec<u8>, String>) -> String {
    names.get(m).cloned().unwrap_or_else(|| {
        let min = forms(m).into_iter().min().unwrap_or_default();
        String::from_utf8_lossy(&min).into_owned()
    })
}

/// `(motif, units)` of the longest tandem run of one unit of `period`
/// bases in `seq`. A run of one base repeated is not a run of a unit.
///
/// After a run the scan restarts one base past where the run's comparison
/// stopped, as the Python does, so a run that starts inside the last partial
/// unit of another can come out up to one unit short.
pub fn longest_run(
    seq: &[u8],
    period: usize,
    names: &HashMap<Vec<u8>, String>,
) -> Option<(String, usize)> {
    let seq = seq.to_ascii_uppercase();
    let n = seq.len();
    let (mut best, mut best_len, mut i) = (None, 0usize, 0usize);
    while i + period < n {
        let mut j = i;
        while j + period < n && seq[j] == seq[j + period] && seq[j] != b'N' {
            j += 1;
        }
        let length = j - i + period;
        let unit = &seq[i..i + period];
        if j > i && length > best_len && unit.iter().any(|&c| c != unit[0]) {
            best = Some(unit);
            best_len = length;
        }
        i = if j > i { j + 1 } else { i + 1 };
    }
    best.map(|u| (motif_name(u, names), best_len / period))
}

const M: u32 = 0;
const INS: u32 = 1;
const DEL: u32 = 2;
const SKIP: u32 = 3;
const SOFT: u32 = 4;
const EQ: u32 = 7;
const DIFF: u32 = 8;

/// The read's bases aligned within `[lo, hi)`, with the soft clip next to
/// them when `softclip` is set; empty when none are.
pub fn repeat_bases(r: &Record, lo: i64, hi: i64, softclip: bool) -> Vec<u8> {
    let seq = r.seq().as_bytes();
    let cig: Vec<(u32, i64)> = r
        .raw_cigar()
        .iter()
        .map(|&c| (c & 0xf, i64::from(c >> 4)))
        .collect();
    let (mut rp, mut q) = (r.pos(), 0i64);
    let (mut first, mut last): (Option<i64>, i64) = (None, 0);
    for &(op, n) in &cig {
        match op {
            M | EQ | DIFF => {
                let (a, b) = (rp.max(lo), (rp + n).min(hi));
                if a < b {
                    first.get_or_insert(q + (a - rp));
                    last = q + (b - rp);
                }
                rp += n;
                q += n;
            }
            DEL | SKIP => rp += n,
            INS | SOFT => {
                if op == INS && first.is_some() && lo <= rp && rp < hi {
                    last = q + n;
                }
                q += n;
            }
            _ => {}
        }
    }
    let Some(mut first) = first else {
        return Vec::new();
    };
    let len = seq.len() as i64;
    if softclip {
        if let Some(&(SOFT, n)) = cig.first()
            && first == n
        {
            first = 0;
        }
        if let Some(&(SOFT, n)) = cig.last()
            && last == len - n
        {
            last = len;
        }
    }
    // Python's slice: clamped to the sequence, empty when reversed.
    let clamp = |x: i64| x.clamp(0, len) as usize;
    let (a, b) = (clamp(first), clamp(last));
    if a < b {
        seq[a..b].to_vec()
    } else {
        Vec::new()
    }
}

fn usable(r: &Record) -> bool {
    r.flags() & (0x400 | 0x100 | 0x800 | 0x200) == 0
}

/// pysam's `AlignmentFile.mate`: the first read at the mate's position with
/// the same name and the other of read 1 and read 2; `None` when there is
/// none. A mate contig past the header has no mate either: htslib answers
/// a query there with no reads, as it does for pysam, whose `ValueError` the
/// Python maps to `None`.
fn mate(bam: &mut Bam, r: &Record) -> Result<Option<Record>, String> {
    let want = (r.flags() ^ 0xC0) & 0xC0;
    if r.mtid() < 0 || r.mpos() < 0 {
        return Ok(None);
    }
    for m in bam
        .fetch_tid(r.mtid(), r.mpos(), r.mpos() + 1)
        .map_err(|e| e.to_string())?
    {
        let m = m.map_err(|e| e.to_string())?;
        if m.flags() & want != 0 && m.qname() == r.qname() {
            return Ok(Some(m));
        }
    }
    Ok(None)
}

/// How the reads of one locus are taken.
struct Window {
    lo: i64,
    hi: i64,
    fetched: (i64, i64),
    softclip: bool,
}

/// Counts the motifs of the reads `fetch` returned, once per read: the
/// bases over the repeat for a mapped read, the whole sequence for an
/// unmapped read and for a mate mapped outside the fetched window
/// (`mates`, `None` turns the rescue off).
fn count_motifs(
    reads: impl Iterator<Item = Result<Record, String>>,
    w: &Window,
    mut mates: Option<&mut Bam>,
    units_min: i64,
    spec: &Spec,
) -> Result<IndexMap<String, usize>, String> {
    let names = spec.names();
    let mut counts: IndexMap<String, usize> = IndexMap::new();
    let mut count = |bases: &[u8]| {
        if let Some((motif, units)) = longest_run(bases, spec.period, &names)
            && units as i64 >= units_min
        {
            *counts.entry(motif).or_insert(0) += 1;
        }
    };
    let mut seen: HashSet<(Vec<u8>, bool)> = HashSet::new();
    for r in reads {
        let r = r?;
        if !usable(&r) {
            continue;
        }
        let read1 = r.is_first_in_template();
        if !seen.insert((r.qname().to_vec(), read1)) {
            continue;
        }
        if r.is_unmapped() {
            count(&r.seq().as_bytes());
            continue;
        }
        count(&repeat_bases(&r, w.lo, w.hi, w.softclip));
        let Some(bam) = mates.as_deref_mut() else {
            continue;
        };
        if !r.is_paired() || r.is_mate_unmapped() {
            continue;
        }
        let mkey = (r.qname().to_vec(), !read1);
        if seen.contains(&mkey) {
            continue;
        }
        if r.mtid() == r.tid() && w.fetched.0 <= r.mpos() && r.mpos() < w.fetched.1 {
            continue; // the mate is in the fetched window: the fetch returns it
        }
        if let Some(m) = mate(bam, &r)?
            && usable(&m)
        {
            seen.insert(mkey);
            count(&m.seq().as_bytes());
        }
    }
    Ok(counts)
}

/// Units a run needs to count: more than the reference allele of `[start, end)`.
fn min_units(start: i64, end: i64, period: usize) -> i64 {
    MIN_UNITS.max((end - start).div_euclid(period as i64) + 1)
}

/// What the motif counts say.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MotifVerdict {
    PathogenicMotif,
    BenignMotif,
    UncertainMotif,
    NotResolved,
    NotChecked,
}

impl MotifVerdict {
    fn name(self) -> &'static str {
        match self {
            MotifVerdict::PathogenicMotif => "pathogenic_motif",
            MotifVerdict::BenignMotif => "benign_motif",
            MotifVerdict::UncertainMotif => "uncertain_motif",
            MotifVerdict::NotResolved => "not_resolved",
            MotifVerdict::NotChecked => "not_checked",
        }
    }
}

/// The columns from `dominant_motif` to `note`.
fn verdict(counts: &IndexMap<String, usize>, units_min: i64, spec: &Spec) -> [String; 6] {
    let total: usize = counts.values().sum();
    let mut order: Vec<(&String, usize)> = counts.iter().map(|(m, &c)| (m, c)).collect();
    order.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    let ratio = |c: usize| c as f64 / total as f64;
    let shares = if total == 0 {
        ".".to_owned()
    } else {
        order
            .iter()
            .map(|(m, c)| format!("{m}={}", fixed(ratio(*c), 2)))
            .collect::<Vec<_>>()
            .join(",")
    };
    let locus = &spec.locus;
    let row = |dom: &str, cls: &str, share: String, v: MotifVerdict, note: String| {
        [
            dom.to_owned(),
            cls.to_owned(),
            share,
            shares.clone(),
            v.name().to_owned(),
            note,
        ]
    };
    if total < MIN_READS {
        return row(
            ".",
            ".",
            ".".into(),
            MotifVerdict::NotResolved,
            format!(
                "{total} read(s) carry {units_min} or more repeat units at {locus}, fewer than {MIN_READS}: \
             the motif is not resolved; confirm the flag with repeat-primed PCR if symptoms fit"
            ),
        );
    }
    let (dom, n) = order[0];
    let cls = spec.known.get(dom).map_or("unknown", |c| c.name());
    let share = ratio(n);
    // Pathogenic motifs most common first; a tie keeps the order seen.
    let mut path: Vec<(&String, usize)> = counts
        .iter()
        .filter(|(m, _)| spec.known.get(*m) == Some(&Class::Pathogenic))
        .map(|(m, &c)| (m, c))
        .collect();
    path.sort_by(|a, b| b.1.cmp(&a.1));
    let npath: usize = path.iter().map(|p| p.1).sum();
    if npath >= MIN_READS && ratio(npath) >= PATHOGENIC_SHARE {
        let names: Vec<&str> = path.iter().map(|p| p.0.as_str()).collect();
        let both = if spec.recessive {
            format!(
                " ({} needs a pathogenic expansion on both copies)",
                spec.disease
            )
        } else {
            String::new()
        };
        return row(
            dom,
            cls,
            fixed(share, 2),
            MotifVerdict::PathogenicMotif,
            format!(
                "{} in {npath} of {total} repeat reads ({}): confirm with repeat-primed PCR{both}",
                names.join("/"),
                percent(ratio(npath), 0)
            ),
        );
    }
    if cls == "benign" && share >= BENIGN_SHARE {
        return row(
            dom,
            cls,
            fixed(share, 2),
            MotifVerdict::BenignMotif,
            format!(
                "motif {dom} (benign) in {n} of {total} repeat reads ({}): the flag is not supported",
                percent(share, 0)
            ),
        );
    }
    row(
        dom,
        cls,
        fixed(share, 2),
        MotifVerdict::UncertainMotif,
        format!(
            "most common motif {dom} ({cls}) in {n} of {total} repeat reads ({}): \
         not a known benign expansion; confirm with repeat-primed PCR if symptoms fit",
            percent(share, 0)
        ),
    )
}

/// A VCF record as the motif count reads it.
#[derive(Debug, Clone)]
struct Rec {
    chrom: String,
    start: i64,
    end: i64,
    status: String,
}

/// REPID and VARID of every record: the first record of a name wins.
fn vcf_records(path: &str) -> Result<HashMap<String, Rec>, String> {
    let f = og_io::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut out: HashMap<String, Rec> = HashMap::new();
    for line in Lines::new(f, Utf8::Strict) {
        let line = line.map_err(|e| format!("{path}: {e}"))?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let mut info: HashMap<&str, &str> = HashMap::new();
        for x in f[7].split(';') {
            let (k, v) = x.split_once('=').unwrap_or((x, ""));
            info.insert(k, v);
        }
        let int = |s: &str| {
            parse_int(s)
                .ok()
                .and_then(|v| i64::try_from(v).ok())
                .ok_or_else(|| format!("{path}: not an integer: {s:?}"))
        };
        let start = int(f[1])? - 1;
        let end = match info.get("END") {
            Some(e) if isdigit(e) => int(e)?,
            _ => start + f[3].chars().count() as i64,
        };
        for key in [info.get("REPID"), info.get("VARID")].into_iter().flatten() {
            if !key.is_empty() && !out.contains_key(*key) {
                out.insert(
                    (*key).to_owned(),
                    Rec {
                        chrom: f[0].to_owned(),
                        start,
                        end,
                        status: info.get("STR_STATUS").copied().unwrap_or("").to_owned(),
                    },
                );
            }
        }
    }
    Ok(out)
}

/// `chr:start-end` of a catalog entry: everything before the first `:`.
fn region(r: &str) -> Option<(String, i64, i64)> {
    let (contig, rest) = r.split_once(':')?;
    let (s, e) = rest.split_once('-')?;
    let num = |x: &str| {
        (!x.is_empty() && x.chars().all(is_decimal))
            .then(|| parse_int(x).ok().and_then(|v| i64::try_from(v).ok()))
            .flatten()
    };
    (!contig.is_empty()).then_some(())?;
    Some((contig.to_owned(), num(s)?, num(e)?))
}

fn strs<'a>(v: &'a Value, what: &str) -> Result<Vec<&'a str>, String> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Str(s) => Ok(s.as_str()),
                _ => Err(format!(
                    "catalog: {what} holds a value that is not a string"
                )),
            })
            .collect(),
        _ => Err(format!("catalog: {what} is not a list")),
    }
}

/// LocusId -> how to read it, for every catalog entry with motifs.
fn catalog_specs(path: &str) -> Result<IndexMap<String, Spec>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let cat = json::loads(&text).map_err(|e| format!("{path}: {e}"))?;
    let Value::Array(entries) = cat else {
        return Err(format!("{path}: not a list of catalog entries"));
    };
    let mut out = IndexMap::new();
    for e in &entries {
        let Value::Object(o) = e else {
            return Err(format!("{path}: an entry is not an object"));
        };
        let Some(path_motifs) = o.get("PathogenicMotifs").filter(|v| truthy(v)) else {
            continue;
        };
        let path_motifs = strs(path_motifs, "PathogenicMotifs")?;
        let mut known: IndexMap<String, Class> = IndexMap::new();
        if let Some(b) = o.get("BenignMotifs") {
            for m in strs(b, "BenignMotifs")? {
                known.insert(m.to_owned(), Class::Benign);
            }
        }
        for &m in &path_motifs {
            known.insert(m.to_owned(), Class::Pathogenic);
        }
        let Some(Value::Str(id)) = o.get("LocusId") else {
            return Err(format!("{path}: an entry has no LocusId string"));
        };
        let regions: Vec<&str> = match o.get("ReferenceRegion") {
            Some(Value::Str(s)) => vec![s.as_str()],
            Some(v @ Value::Array(_)) => strs(v, "ReferenceRegion")?,
            _ => return Err(format!("{path}: {id} has no ReferenceRegion")),
        };
        let regions = regions
            .iter()
            .map(|r| region(r).ok_or_else(|| format!("{path}: {id}: not a region: {r:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        let Some(first) = regions.first() else {
            return Err(format!("{path}: {id} has no ReferenceRegion"));
        };
        let span = (
            first.0.clone(),
            regions.iter().map(|r| r.1).min().unwrap_or(0),
            regions.iter().map(|r| r.2).max().unwrap_or(0),
        );
        let disease = match o.get("Disease") {
            None => id.clone(),
            Some(v) => py_str(v).map_err(|e| format!("{path}: {id}: Disease: {e:?}"))?,
        };
        let recessive = match o.get("InheritanceMode") {
            None => false,
            Some(Value::Str(s)) => s == "AR",
            Some(_) => return Err(format!("{path}: {id}: InheritanceMode is not a string")),
        };
        out.insert(
            id.clone(),
            Spec {
                locus: id.clone(),
                known,
                period: path_motifs[0].chars().count(),
                disease,
                recessive,
                span: Some(span),
            },
        );
    }
    Ok(out)
}

fn flagged(status: &str) -> bool {
    status
        .split(',')
        .any(|x| x == "pre_mutation" || x == "full_mutation")
}

/// What is asked of one run.
struct Opts<'a> {
    bam: Option<&'a str>,
    skip: Option<&'a str>,
    softclip: bool,
    mates: bool,
}

/// The output row of one locus.
fn check(o: &Opts<'_>, spec: &Spec, records: &HashMap<String, Rec>) -> Result<Vec<String>, String> {
    let name = &spec.locus;
    let Some(rec) = records.get(name) else {
        return Ok([
            name.as_str(),
            ".",
            ".",
            ".",
            "0",
            ".",
            ".",
            ".",
            ".",
            "not_checked",
        ]
        .iter()
        .map(|s| (*s).to_owned())
        .chain([format!("no {name} record in the VCF")])
        .collect());
    };
    let (chrom, start, end) = spec
        .span
        .clone()
        .unwrap_or_else(|| (rec.chrom.clone(), rec.start, rec.end));
    let mut row = vec![
        name.clone(),
        chrom.clone(),
        (start + 1).to_string(),
        end.to_string(),
    ];
    let empty = |row: &mut Vec<String>, note: String| {
        row.extend(["0", ".", ".", ".", ".", "not_checked"].map(String::from));
        row.push(note);
    };
    if let Some(reason) = o.skip {
        empty(&mut row, reason.to_owned());
        return Ok(row);
    }
    let bam_path = o.bam.ok_or("give --bam, or --skip REASON")?;
    let open = || Bam::open(Path::new(bam_path)).map_err(|e| format!("{bam_path}: {e}"));
    let mut bam = open()?;
    if bam.tid(&chrom).is_none() {
        empty(&mut row, format!("the BAM has no contig {chrom}"));
        return Ok(row);
    }
    let mut other = if o.mates { Some(open()?) } else { None };
    let w = Window {
        lo: start - FLANK,
        hi: end + FLANK,
        fetched: ((start - FETCH_PAD).max(0), end + FETCH_PAD),
        softclip: o.softclip,
    };
    let units = min_units(start, end, spec.period);
    let reads = bam
        .fetch(&chrom, w.fetched.0, w.fetched.1)
        .map_err(|e| format!("{bam_path}: {e}"))?
        .map(|r| r.map_err(|e| format!("{bam_path}: {e}")));
    let counts = count_motifs(reads, &w, other.as_mut(), units, spec)?;
    row.push(counts.values().sum::<usize>().to_string());
    row.extend(verdict(&counts, units, spec));
    Ok(row)
}

fn run(a: &Args) -> Result<String, String> {
    let (vcf, out) = (a.required("--vcf")?, a.required("--out")?);
    let o = Opts {
        bam: a.get("--bam"),
        // `--skip ''` reads as no reason, as in the Python.
        skip: a.get("--skip").filter(|s| !s.is_empty()),
        softclip: !a.flag("--no-softclip"),
        mates: !a.flag("--no-mates"),
    };
    let specs = match a.get("--catalog") {
        Some(c) => catalog_specs(c)?,
        None => IndexMap::new(),
    };
    let records = vcf_records(vcf)?;
    let todo: Vec<Spec> = if a.flag("--flagged") {
        specs
            .into_iter()
            .filter(|(n, _)| n != "RFC1" && records.get(n).is_some_and(|r| flagged(&r.status)))
            .map(|(_, s)| s)
            .collect()
    } else {
        vec![specs.get("RFC1").cloned().unwrap_or_else(Spec::rfc1)]
    };
    let rows = todo
        .iter()
        .map(|s| check(&o, s, &records))
        .collect::<Result<Vec<_>, _>>()?;
    if rows.is_empty() {
        match std::fs::remove_file(out) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                return Err(format!("cannot remove {out}: {e}"));
            }
            _ => {}
        }
        return Ok("No flagged locus with motifs in the catalog: no motif count\n".to_owned());
    }
    let mut text = format!("{COLUMNS}\n");
    let mut said = String::new();
    for r in &rows {
        text.push_str(&r.join("\t"));
        text.push('\n');
        said.push_str(&format!("{} motif count: {} ({})\n", r[0], r[9], r[10]));
    }
    let mut f = AtomicFile::create(out).map_err(|e| format!("cannot write {out}: {e}"))?;
    f.write_all(text.as_bytes())
        .and_then(|()| f.commit())
        .map_err(|e| format!("cannot write {out}: {e}"))?;
    Ok(said)
}

/// `og str motif ARGS`: the exit code; a refused command line is exit 2.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    const USAGE: &str = "usage: og str motif --vcf VCF --out OUT [--bam BAM] [--skip SKIP] \
                         [--no-softclip] [--no-mates] [--catalog CATALOG] [--flagged]";
    let parsed = Args::parse(
        argv,
        &["--vcf", "--out", "--bam", "--skip", "--catalog"],
        &["--no-softclip", "--no-mates", "--flagged"],
    )
    .and_then(|a| {
        a.required("--vcf")?;
        a.required("--out")?;
        if a.get("--skip").is_none_or(str::is_empty) && a.get("--bam").is_none() {
            return Err("give --bam, or --skip REASON".to_owned());
        }
        if a.flag("--flagged") && a.get("--catalog").is_none() {
            return Err("--flagged needs --catalog".to_owned());
        }
        Ok(a)
    });
    let a = match parsed {
        Ok(a) => a,
        Err(e) => {
            let _ = writeln!(err, "{USAGE}\nog: error: {e}");
            return 2;
        }
    };
    match run(&a).and_then(|s| out.write_all(s.as_bytes()).map_err(|e| e.to_string())) {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rfc1_names() -> HashMap<Vec<u8>, String> {
        Spec::rfc1().names()
    }

    #[test]
    fn a_motif_is_named_by_its_known_form_on_either_strand() {
        let n = rfc1_names();
        assert_eq!(motif_name(b"GGGAA", &n), "AAGGG");
        assert_eq!(motif_name(b"CCCTT", &n), "AAGGG");
        assert_eq!(motif_name(b"GAGAA", &n), "AAGAG");
    }

    #[test]
    fn longest_run_finds_the_unit_and_skips_homopolymers() {
        let n = rfc1_names();
        let seq = [b"ACGTTGCA".as_slice(), &b"AAAAG".repeat(11), b"TTGCACGT"].concat();
        assert_eq!(longest_run(&seq, 5, &n), Some(("AAAAG".into(), 11)));
        assert_eq!(longest_run(&b"A".repeat(60), 5, &n), None);
        // Reverse strand: CCCTT is AAGGG.
        assert_eq!(
            longest_run(&b"CCCTT".repeat(8), 5, &n).map(|r| r.0),
            Some("AAGGG".into())
        );
    }

    #[test]
    fn forms_has_every_rotation_on_both_strands() {
        let f = forms(b"AAGGG");
        assert_eq!(f.len(), 10);
        assert!(f.contains(&b"CCCTT".to_vec()) && f.contains(&b"GGAAG".to_vec()));
    }

    #[test]
    fn the_verdicts_and_their_words() {
        let s = Spec::rfc1();
        let c = |v: &[(&str, usize)]| -> IndexMap<String, usize> {
            v.iter().map(|(m, n)| ((*m).to_owned(), *n)).collect()
        };
        let v = verdict(&c(&[("AAGGG", 13), ("AAAAG", 3)]), 12, &s);
        assert_eq!(v[3], "AAGGG=0.81,AAAAG=0.19");
        assert_eq!(v[4], "pathogenic_motif");
        assert_eq!(
            v[5],
            "AAGGG in 13 of 16 repeat reads (81%): confirm with repeat-primed PCR \
             (CANVAS needs a pathogenic expansion on both copies)"
        );
        assert_eq!(
            verdict(&c(&[("AAAAG", 6), ("AAAGG", 4)]), 12, &s)[4],
            "benign_motif"
        );
        assert_eq!(
            verdict(&c(&[("AAGAG", 5)]), 12, &s)[..2],
            ["AAGAG", "unknown"]
        );
        assert_eq!(verdict(&c(&[("AAGGG", 2)]), 12, &s)[4], "not_resolved");
        assert_eq!(verdict(&c(&[]), 12, &s)[3], ".");
    }

    #[test]
    fn the_floor_is_one_unit_more_than_the_reference_allele() {
        assert_eq!(min_units(39_348_424, 39_348_479, 5), 12);
        assert_eq!(min_units(0, 5, 5), 4);
    }

    #[test]
    fn a_mate_on_a_contig_the_header_lacks_is_no_mate() {
        use rust_htslib::bam::{self as hb, header::HeaderRecord, record::Cigar};
        let dir = std::env::temp_dir().join(format!("og-bam-mate-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut h = hb::Header::new();
        h.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr4")
                .push_tag(b"LN", 1000),
        );
        let write = |name: &str, mtid: i32| {
            let path = dir.join(name);
            let mut rec = Record::new();
            let cigar = hb::record::CigarString(vec![Cigar::Match(4)]);
            rec.set(b"r1", Some(&cigar), b"ACGT", &[30; 4]);
            rec.set_tid(0);
            rec.set_pos(10);
            rec.set_mtid(mtid);
            rec.set_mpos(10);
            rec.set_flags(0x1 | 0x40);
            let mut w = hb::Writer::from_path(&path, &h, hb::Format::Bam).unwrap();
            w.set_compression_level(hb::CompressionLevel::Uncompressed)
                .unwrap();
            w.write(&rec).unwrap();
            path
        };
        // htslib will not index a record whose mate contig is past the
        // header, so the index comes from the same bytes with mtid 0.
        let good = write("good.bam", 0);
        let bad = write("bad.bam", 5);
        hb::index::build(&good, None, hb::index::Type::Bai, 1).unwrap();
        std::fs::copy(dir.join("good.bam.bai"), dir.join("bad.bam.bai")).unwrap();
        let mut bam = Bam::open(&bad).unwrap();
        let read = bam.fetch("chr4", 0, 100).unwrap().next().unwrap().unwrap();
        assert_eq!(read.mtid(), 5);
        assert!(matches!(mate(&mut bam, &read), Ok(None)));
        drop(bam);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
