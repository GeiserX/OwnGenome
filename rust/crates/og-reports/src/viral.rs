//! `og viral report` (step 47): the viral DNA in the reads that do not map
//! to the human genome.
//!
//! ```text
//! og viral report --sample S --viral-fasta viral_genomes.fa
//!                 --viral-coverage S_viral_coverage.tsv --viral-depth S_viral_depth.tsv
//!                 --human-coverage S_human_coverage.tsv --viral-hits S_viral_hits.tsv
//!                 --idxstats S_idxstats.tsv --out S_viral.tsv
//! ```
//!
//! The unmapped reads are realigned to the viral genomes; `samtools
//! coverage` of that BAM gives each genome's length and aligned records,
//! `samtools depth -aa` the depth at every viral base, `samtools view -F 4`
//! of it the alignments themselves, and `samtools coverage` of the sample's
//! BAM the human depth (chr1 to chr22, over the bases with at least one
//! read). `samtools idxstats` of the sample's BAM says whether it still holds
//! its unplaced unmapped reads. One row per viral genome, in the FASTA's order.
//!
//! The evidence for a call is the *counted* reads: primary alignments with a
//! mapping quality of 20 or more, not low complexity (a DUST score above 7,
//! see [`dust`]) and not a duplicate (the same genome, strand and unclipped
//! 5' end as a read already counted). `reads` and `breadth` are theirs;
//! supplementary records of counted quality add breadth, not reads. The
//! integration rule reads the depth of every primary alignment instead.
//!
//! The call: `not assessed` for every genome when the BAM has no unplaced
//! unmapped read (cut to regions or filtered, so nothing could be found);
//! inherited integrated HHV-6 when the mode of the per-base depth over
//! HHV-6A or HHV-6B is at least 0.30 of the human depth (Liu et al. 2020,
//! PLoS Genet, doi:10.1371/journal.pgen.1008915: one integrated copy per
//! cell gives about 0.5, viraemia 10 to 1000 times less; when both species
//! pass, both are called); `lab control (phiX)` for Illumina's spike-in;
//! `viral DNA found` (10 counted reads or more over 10% of the genome or
//! more); `a few reads` (10 or more over less); `below threshold` (1 to 9
//! counted reads); `not found` (none). None of these is a diagnosis of
//! infection.
//!
//! The depth file (about 1.6 million lines over the 24 genomes) is read
//! once into a histogram per genome: memory is bounded by the number of
//! distinct depths, not the genome length. The alignments are read once,
//! one line at a time, into a bit per genome base and a set of the read
//! starts already counted.
//!
//! What changes on purpose: a number in a samtools table that is not one
//! stops with the file and line (`ERROR:`, exit 1) instead of a traceback,
//! and problems start `ERROR:`, not `viral_report.py:`.

use crate::Error;
use indexmap::IndexMap;
use og_filters::args::Spec;
use og_io::Utf8;
use rustc_hash::FxHashSet;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

/// Counted reads a call needs; fewer read `below threshold`.
const MIN_READS: i128 = 10;
/// `viral DNA found` needs the counted reads over a tenth of the genome:
/// `covbases * FOUND_BREADTH_DIV >= length`, in whole numbers.
const FOUND_BREADTH_DIV: i128 = 10;
const ICI_RELATIVE_DEPTH: f64 = 0.30;
/// The lowest mapping quality a counted read has.
const MIN_MAPQ: u32 = 20;
/// A read whose DUST score is above this is low complexity and not counted.
const DUST_MAX: f64 = 7.0;
pub const HEADER: [&str; 12] = [
    "virus",
    "accession",
    "length",
    "reads",
    "breadth",
    "mean_depth",
    "mode_depth",
    "human_depth",
    "relative_depth",
    "copies_per_cell",
    "call",
    "note",
];

/// The call of one viral genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViralCall {
    Ici,
    LabControl,
    Found,
    FewReads,
    BelowThreshold,
    NotFound,
    NotAssessed,
}

impl ViralCall {
    fn text(self) -> &'static str {
        match self {
            ViralCall::Ici => "inherited integrated HHV-6 likely",
            ViralCall::LabControl => "lab control (phiX)",
            ViralCall::Found => "viral DNA found",
            ViralCall::FewReads => "a few reads",
            ViralCall::BelowThreshold => "below threshold",
            ViralCall::NotFound => "not found",
            ViralCall::NotAssessed => "not assessed",
        }
    }
}

fn hhv6(base: &str) -> Option<&'static str> {
    match base {
        "NC_001664" => Some("HHV-6A"),
        "NC_000898" => Some("HHV-6B"),
        _ => None,
    }
}

const PHIX: &str = "NC_001422";
/// HCV: an RNA virus with no DNA stage.
const RNA_VIRUS: &str = "NC_004102";

fn base(acc: &str) -> &str {
    acc.split('.').next().unwrap_or(acc)
}

/// `[(accession, title)]` of the FASTA, in file order.
fn read_fasta_titles(path: &str) -> Result<Vec<(String, String)>, Error> {
    let mut out = Vec::new();
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(file, Utf8::Replace);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        if let Some(rest) = line.strip_prefix('>') {
            let t = crate::strip(rest);
            let (acc, title) = t.split_once(' ').unwrap_or((t, ""));
            let title = crate::strip(title);
            out.push((
                acc.to_owned(),
                if title.is_empty() { acc } else { title }.to_owned(),
            ));
        }
    }
    if out.is_empty() {
        return Err(Error(format!("{path} has no sequence")));
    }
    Ok(out)
}

/// One `samtools coverage` row.
#[derive(Debug, Clone, Copy)]
struct Cov {
    reads: i128,
    length: i128,
    covbases: i128,
    meandepth: f64,
}

/// `samtools coverage` rows by contig, in file order.
fn read_coverage(path: &str) -> Result<IndexMap<String, Cov>, Error> {
    let mut rows = IndexMap::new();
    let mut header: Option<Vec<String>> = None;
    for (i, line) in crate::read_lines(path, Utf8::Replace)?.iter().enumerate() {
        let cols: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(line)
            .split('\t')
            .collect();
        if line.starts_with('#') {
            header = Some(
                cols.iter()
                    .map(|c| c.trim_start_matches('#').to_owned())
                    .collect(),
            );
            continue;
        }
        let Some(h) = header.as_ref().filter(|h| h.len() == cols.len()) else {
            return Err(Error(format!("{path} is not samtools coverage output")));
        };
        let r: HashMap<&str, &str> = h
            .iter()
            .map(String::as_str)
            .zip(cols.iter().copied())
            .collect();
        let cell = |k: &str| {
            r.get(k).copied().ok_or_else(|| {
                Error(format!(
                    "{path} is not samtools coverage output (no {k} column)"
                ))
            })
        };
        let bad = |k: &str, v: &str| {
            Error(format!(
                "{path} line {}: {k} {} is not a number",
                i + 1,
                og_filters::py_repr(v)
            ))
        };
        let int = |k: &str| -> Result<i128, Error> {
            let v = cell(k)?;
            crate::int(v).ok_or_else(|| bad(k, v))
        };
        let mean = cell("meandepth")?;
        let cov = Cov {
            reads: int("numreads")?,
            length: int("endpos")? - int("startpos")? + 1,
            covbases: int("covbases")?,
            meandepth: crate::float(mean).ok_or_else(|| bad("meandepth", mean))?,
        };
        rows.insert(cell("rname")?.to_owned(), cov);
    }
    Ok(rows)
}

/// chr1 to chr22, with or without the prefix.
fn autosome(name: &str) -> bool {
    let n = name.strip_prefix("chr").unwrap_or(name);
    let b = n.as_bytes();
    match b {
        [d] => (b'1'..=b'9').contains(d),
        [b'1', d] => d.is_ascii_digit(),
        [b'2', d] => (b'0'..=b'2').contains(d),
        _ => false,
    }
}

/// Mean depth over chr1 to chr22's bases with at least one read.
fn human_depth(cov: &IndexMap<String, Cov>) -> Option<f64> {
    let (mut total, mut covered) = (0.0f64, 0.0f64);
    for (name, r) in cov {
        if autosome(name) && r.covbases > 0 {
            total += r.meandepth * r.length as f64;
            covered += r.covbases as f64;
        }
    }
    (covered != 0.0).then(|| total / covered)
}

/// Depths seen over one genome: small depths in a vector, the rest in a map.
#[derive(Default)]
struct Hist {
    small: Vec<u64>,
    big: BTreeMap<i128, u64>,
    positions: u64,
}

const SMALL: usize = 4096;

impl Hist {
    fn add(&mut self, d: i128) {
        self.positions += 1;
        match usize::try_from(d) {
            Ok(i) if i < SMALL => {
                if self.small.len() <= i {
                    self.small.resize(i + 1, 0);
                }
                self.small[i] += 1;
            }
            _ => *self.big.entry(d).or_default() += 1,
        }
    }

    /// (depth, count) of every depth seen, smallest first.
    fn counts(&self) -> Vec<(i128, u64)> {
        let mut v: Vec<(i128, u64)> = self.big.range(..0).map(|(d, c)| (*d, *c)).collect();
        v.extend(
            self.small
                .iter()
                .enumerate()
                .filter(|(_, c)| **c > 0)
                .map(|(d, c)| (d as i128, *c)),
        );
        v.extend(self.big.range(0..).map(|(d, c)| (*d, *c)));
        v
    }
}

/// (mean, mode) per genome of `lengths`, from `samtools depth -aa`.
fn read_depth(
    path: &str,
    lengths: &IndexMap<String, i128>,
) -> Result<HashMap<String, (f64, i128)>, Error> {
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(file, Utf8::Replace);
    let mut hist: HashMap<String, Hist> = HashMap::new();
    let mut line = String::new();
    let mut n = 0u64;
    let mut last: Option<String> = None;
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        n += 1;
        let mut cols = line.split('\t');
        let (Some(name), Some(_), Some(d), None) =
            (cols.next(), cols.next(), cols.next(), cols.next())
        else {
            return Err(Error(format!(
                "{path} line {n} is not samtools depth output (contig, position, depth)"
            )));
        };
        let d = crate::int(d).ok_or_else(|| {
            Error(format!(
                "{path} line {n}: depth {} is not a whole number",
                og_filters::py_repr(d.trim_end_matches('\n'))
            ))
        })?;
        // The file runs contig by contig: look the contig up only when it changes.
        if last.as_deref() != Some(name) {
            last = Some(name.to_owned());
            hist.entry(name.to_owned()).or_default();
        }
        if let Some(h) = hist.get_mut(name) {
            h.add(d);
        }
    }
    let mut out = HashMap::new();
    for (name, &length) in lengths {
        let seen = hist.get(name).map_or(0, |h| h.positions);
        if seen != 0 && i128::from(seen) != length {
            return Err(Error(format!(
                "{path} has {seen} positions for {name}, its length is {length} (run samtools depth -aa)"
            )));
        }
        if length == 0 {
            return Err(Error(format!("{name} has length 0 in the viral coverage")));
        }
        let counts = match hist.get(name) {
            Some(h) if h.positions > 0 => h.counts(),
            _ => vec![(0, u64::try_from(length).unwrap_or(0))],
        };
        let sum: i128 = counts.iter().map(|(d, c)| d * i128::from(*c)).sum();
        let mean = sum as f64 / length as f64;
        // The most common depth; a tie goes to the lower one.
        let mode = counts
            .iter()
            .fold(None::<(i128, u64)>, |best, &(d, c)| match best {
                Some((_, bc)) if bc >= c => best,
                _ => Some((d, c)),
            })
            .map_or(0, |(d, _)| d);
        out.insert(name.clone(), (mean, mode));
    }
    Ok(out)
}

/// The DUST score of a read, 0 to 100, as PRINSEQ computes it (Schmieder
/// and Edwards 2011, doi:10.1093/bioinformatics/btr026): over windows of 64
/// bases every 32 bases (one more at the read's end when the steps leave a
/// tail; one window for a read under 64 bases), the DUST score of Morgulis
/// et al. 2006 (doi:10.1089/cmb.2006.13.1028) of each window's triplets,
/// `sum c(c-1)/2 / (l-1)`, then their mean scaled by 100/31 (31 is a
/// 64-base homopolymer's score). Random sequence scores about 1.5; a
/// dinucleotide repeat about 49, a homopolymer 100. Triplets with a base
/// other than A, C, G or T are left out.
pub fn dust(seq: &[u8]) -> f64 {
    const WINDOW: usize = 64;
    const STEP: usize = 32;
    let n = seq.len();
    let mut starts: Vec<usize> = (0..)
        .step_by(STEP)
        .take_while(|s| s + WINDOW <= n)
        .collect();
    match starts.last() {
        None => starts.push(0),
        Some(&last) if last + WINDOW < n => starts.push(n - WINDOW),
        Some(_) => {}
    }
    let total: f64 = starts
        .iter()
        .map(|&s| window_dust(&seq[s..(s + WINDOW).min(n)]))
        .sum();
    total / starts.len() as f64 * 100.0 / 31.0
}

/// The DUST score of one window: `sum c(c-1)/2 / (l-1)` over its triplets.
fn window_dust(w: &[u8]) -> f64 {
    let code = |b: u8| match b.to_ascii_uppercase() {
        b'A' => Some(0usize),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    };
    let mut counts = [0u32; 64];
    let mut l = 0u32;
    for t in w.windows(3) {
        if let (Some(a), Some(b), Some(c)) = (code(t[0]), code(t[1]), code(t[2])) {
            counts[a * 16 + b * 4 + c] += 1;
            l += 1;
        }
    }
    if l < 2 {
        return 0.0;
    }
    let pairs: u32 = counts.iter().map(|&c| c * c.saturating_sub(1) / 2).sum();
    f64::from(pairs) / f64::from(l - 1)
}

/// What the counted reads of one genome give.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Counted {
    /// Counted primary reads.
    reads: i128,
    /// Bases covered by a counted record (primary or supplementary).
    covbases: i128,
}

/// One alignment line of `samtools view`: `None` for a record that never
/// counts, else (contig, reverse, unclipped 5' end, supplementary, blocks).
type Hit<'a> = (&'a str, bool, i64, bool, Vec<(i64, i64)>);

/// Reads one SAM line; `Err` names what is wrong with it.
fn parse_hit(line: &str) -> Result<Option<Hit<'_>>, String> {
    let f: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
    if f.len() < 11 {
        return Err("is not a SAM record (11 columns or more)".into());
    }
    let flag: u16 = f[1]
        .parse()
        .map_err(|_| format!("FLAG {} is not a number", f[1]))?;
    // Unmapped, secondary, QC-failed or marked duplicate: never counted.
    if flag & (0x4 | 0x100 | 0x200 | 0x400) != 0 {
        return Ok(None);
    }
    let mapq: u32 = f[4]
        .parse()
        .map_err(|_| format!("MAPQ {} is not a number", f[4]))?;
    if mapq < MIN_MAPQ || mapq == 255 {
        return Ok(None);
    }
    if f[9] != "*" && dust(f[9].as_bytes()) > DUST_MAX {
        return Ok(None);
    }
    let pos: i64 = f[3]
        .parse()
        .map_err(|_| format!("POS {} is not a number", f[3]))?;
    let mut at = pos;
    let (mut lead, mut trail, mut blocks, mut seen_ref) = (0i64, 0i64, Vec::new(), false);
    let mut num = String::new();
    for c in f[5].chars() {
        if c.is_ascii_digit() {
            num.push(c);
            continue;
        }
        let n: i64 = num
            .parse()
            .map_err(|_| format!("CIGAR {} is not a CIGAR", f[5]))?;
        num.clear();
        match c {
            'M' | '=' | 'X' => {
                blocks.push((at, at + n));
                at += n;
                seen_ref = true;
                trail = 0;
            }
            'D' | 'N' => {
                at += n;
                seen_ref = true;
                trail = 0;
            }
            'S' | 'H' => {
                if seen_ref {
                    trail += n;
                } else {
                    lead += n;
                }
            }
            'I' | 'P' => {}
            _ => return Err(format!("CIGAR {} is not a CIGAR", f[5])),
        }
    }
    if !num.is_empty() || blocks.is_empty() {
        return Err(format!("CIGAR {} has no aligned base", f[5]));
    }
    let reverse = flag & 0x10 != 0;
    let five = if reverse { at - 1 + trail } else { pos - lead };
    Ok(Some((f[2], reverse, five, flag & 0x800 != 0, blocks)))
}

/// The counted reads per genome of `lengths`, from `samtools view -F 4` of
/// the viral BAM (header lines are skipped).
fn read_hits(
    path: &str,
    lengths: &IndexMap<String, i128>,
) -> Result<HashMap<String, Counted>, Error> {
    let mut covered: Vec<Vec<u64>> = vec![Vec::new(); lengths.len()];
    let mut reads = vec![0i128; lengths.len()];
    let mut seen: FxHashSet<(usize, bool, i64, bool)> = FxHashSet::default();
    for (n, line) in crate::lines(path, Utf8::Replace)?.enumerate() {
        let line = line?;
        if line.starts_with('@') || line.trim().is_empty() {
            continue;
        }
        let hit = parse_hit(&line).map_err(|e| Error(format!("{path} line {}: {e}", n + 1)))?;
        let Some((contig, reverse, five, supplementary, blocks)) = hit else {
            continue;
        };
        let Some((g, _, &length)) = lengths.get_full(contig) else {
            return Err(Error(format!(
                "{path} line {}: {contig} is not in the viral coverage: was the BAM aligned to this FASTA?",
                n + 1
            )));
        };
        if !seen.insert((g, reverse, five, supplementary)) {
            continue;
        }
        if !supplementary {
            reads[g] += 1;
        }
        let bits = &mut covered[g];
        if bits.is_empty() {
            bits.resize(usize::try_from(length).unwrap_or(0).div_ceil(64), 0);
        }
        for (from, to) in blocks {
            for p in from.max(1)..to.min(i64::try_from(length).unwrap_or(i64::MAX) + 1) {
                let i = usize::try_from(p - 1).unwrap_or(0);
                bits[i / 64] |= 1 << (i % 64);
            }
        }
    }
    Ok(lengths
        .keys()
        .enumerate()
        .map(|(g, name)| {
            let covbases = covered[g].iter().map(|w| i128::from(w.count_ones())).sum();
            (
                name.clone(),
                Counted {
                    reads: reads[g],
                    covbases,
                },
            )
        })
        .collect())
}

/// The unplaced unmapped reads of the sample's BAM: the fourth column of
/// the `*` row of `samtools idxstats`.
fn read_unplaced(path: &str) -> Result<i128, Error> {
    let mut found = None;
    for (n, line) in crate::lines(path, Utf8::Replace)?.enumerate() {
        let line = line?;
        let cols: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
        if cols.len() != 4 {
            if line.trim().is_empty() {
                continue;
            }
            return Err(Error(format!(
                "{path} line {} is not samtools idxstats output (contig, length, mapped, unmapped)",
                n + 1
            )));
        }
        if cols[0] == "*" {
            found = Some(crate::int(cols[3]).ok_or_else(|| {
                Error(format!(
                    "{path} line {}: unmapped {} is not a number",
                    n + 1,
                    og_filters::py_repr(cols[3])
                ))
            })?);
        }
    }
    found.ok_or_else(|| {
        Error(format!(
            "{path} has no * row: is it samtools idxstats output?"
        ))
    })
}

/// A word character for Python's `\b` on a str.
fn word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// The title without its `, complete genome` tail:
/// `re.sub(r",? (complete (virion )?genome|genome|DNA)\b.*$", "", title)`.
fn short_title(title: &str) -> &str {
    const ALTS: [&str; 4] = ["complete virion genome", "complete genome", "genome", "DNA"];
    let at = |rest: &str| {
        ALTS.iter().any(|a| {
            rest.strip_prefix(a)
                .is_some_and(|after| !after.chars().next().is_some_and(word))
        })
    };
    for (i, _) in title.char_indices() {
        let rest = &title[i..];
        let tail = rest.strip_prefix(", ").or_else(|| rest.strip_prefix(' '));
        if tail.is_some_and(at) {
            return &title[..i];
        }
    }
    title
}

fn short_name(acc: &str, title: &str) -> String {
    if let Some(n) = hhv6(base(acc)) {
        return n.into();
    }
    let t = crate::strip(short_title(title));
    if t.is_empty() { acc } else { t }.to_owned()
}

/// What one genome's call reads.
#[derive(Debug, Clone, Copy)]
struct Evidence {
    /// Every aligned record `samtools coverage` counts.
    aligned: i128,
    /// The counted reads and the bases they cover.
    counted: Counted,
    length: i128,
    /// Mode depth over the human depth; `None` without a human depth.
    relative: Option<f64>,
}

/// The call and its note for one viral genome (before the HHV-6 pair rule).
fn call_row(acc: &str, e: &Evidence) -> (ViralCall, String) {
    let b = base(acc);
    let reads = e.counted.reads;
    if b == PHIX {
        if reads != 0 {
            return (
                ViralCall::LabControl,
                "Illumina adds phiX to the run as a control: these reads were added in the lab."
                    .into(),
            );
        }
        return (
            ViralCall::NotFound,
            "Illumina's phiX control: none of its reads reached this sample's data.".into(),
        );
    }
    if let (Some(_), Some(rel)) = (hhv6(b), e.relative)
        && rel >= ICI_RELATIVE_DEPTH
    {
        return (
            ViralCall::Ici,
            format!(
                "Depth over the genome at {} of the human depth (Liu et al. 2020 call integration at {} or more; one copy per cell gives about 0.5). Inherited integrated HHV-6 is in every cell and is passed on to about half of the children; it can make a blood test read as an active HHV-6 infection.",
                crate::fixed(rel, 2),
                crate::fixed(ICI_RELATIVE_DEPTH, 2)
            ),
        );
    }
    if reads == 0 {
        let note = if e.aligned > 0 {
            format!(
                "{} aligned records, none counted: duplicates, mapping quality under {MIN_MAPQ}, low complexity or supplementary.",
                e.aligned
            )
        } else {
            String::new()
        };
        return (ViralCall::NotFound, note);
    }
    if reads < MIN_READS {
        return (
            ViralCall::BelowThreshold,
            format!(
                "{reads} counted read{}, under the {MIN_READS} a call needs: too few to tell the virus from contamination at this depth.",
                if reads == 1 { "" } else { "s" }
            ),
        );
    }
    let note = if hhv6(b).is_some() {
        if e.relative.is_some() {
            "Below the integration threshold of Liu et al. 2020: not the one-copy-per-cell pattern of inherited HHV-6."
        } else {
            "No human depth (no read on chr1 to chr22): integration is not assessed."
        }
    } else if b == RNA_VIRUS {
        "An RNA virus with no DNA stage: DNA reads point to a contamination or an integration (Moustafa et al. 2017)."
    } else {
        ""
    };
    let call = if e.counted.covbases * FOUND_BREADTH_DIV >= e.length {
        ViralCall::Found
    } else {
        ViralCall::FewReads
    };
    (call, note.into())
}

/// One row: the cells of [`HEADER`] before `copies_per_cell`, and the call.
struct Row {
    virus: String,
    cells: Vec<String>,
    call: ViralCall,
    note: String,
    reads: i128,
    mean: f64,
    mode: i128,
}

/// The files of `og viral report`.
pub struct Files<'a> {
    pub viral_fasta: &'a str,
    pub viral_coverage: &'a str,
    pub viral_depth: &'a str,
    pub human_coverage: &'a str,
    /// `samtools view -F 4` of the viral BAM.
    pub viral_hits: &'a str,
    /// `samtools idxstats` of the sample's BAM.
    pub idxstats: &'a str,
    pub out: &'a str,
}

/// `og viral report`: writes the table; the stdout line.
pub fn report(sample: &str, f: &Files<'_>) -> Result<String, Error> {
    let titles = read_fasta_titles(f.viral_fasta)?;
    let vcov = read_coverage(f.viral_coverage)?;
    let lengths: IndexMap<String, i128> = titles
        .iter()
        .filter_map(|(a, _)| vcov.get(a).map(|c| (a.clone(), c.length)))
        .collect();
    let vdepth = read_depth(f.viral_depth, &lengths)?;
    let hdepth = human_depth(&read_coverage(f.human_coverage)?).filter(|h| *h != 0.0);
    let hits = read_hits(f.viral_hits, &lengths)?;
    let assessed = read_unplaced(f.idxstats)? > 0;
    let mut rows = Vec::new();
    for (acc, title) in &titles {
        let Some(c) = vcov.get(acc) else {
            return Err(Error(format!(
                "the viral coverage has no row for {acc}: was the BAM aligned to this FASTA?"
            )));
        };
        let (mean, mode) = vdepth[acc];
        let counted = hits[acc];
        let relative = hdepth.map(|h| mode as f64 / h);
        let e = Evidence {
            aligned: c.reads,
            counted,
            length: c.length,
            relative,
        };
        let (call, note) = if assessed {
            call_row(acc, &e)
        } else {
            (
                ViralCall::NotAssessed,
                "The BAM holds no unplaced unmapped read (cut to regions or filtered): nothing could be found, so nothing is ruled out.".into(),
            )
        };
        let virus = short_name(acc, title);
        rows.push(Row {
            cells: vec![
                virus.clone(),
                acc.clone(),
                c.length.to_string(),
                counted.reads.to_string(),
                crate::fixed(counted.covbases as f64 / c.length as f64, 4),
                crate::fixed(mean, 3),
                mode.to_string(),
                hdepth.map_or(".".into(), |h| crate::fixed(h, 2)),
                relative.map_or(".".into(), |r| crate::fixed(r, 3)),
            ],
            virus,
            call,
            note,
            reads: counted.reads,
            mean,
            mode,
        });
    }
    // Both HHV-6 species pass: one integrated genome does not lift the
    // other species to a third of the human depth, so both are called.
    let ici: Vec<usize> = (0..rows.len())
        .filter(|&i| rows[i].call == ViralCall::Ici)
        .collect();
    if ici.len() > 1 {
        for &i in &ici {
            rows[i].note.push_str(
                " HHV-6A and HHV-6B both pass: one integrated genome does not bring the other species to this depth, so both are likely integrated; a laboratory test tells them apart.",
            );
        }
    }
    let mut w = pycompat::csv::Writer::new(Vec::new(), '\t');
    let io = |e: std::io::Error| Error::io(&e, f.out);
    w.write_row(HEADER).map_err(io)?;
    for r in &rows {
        // Copies per cell where the call stands on depth: the mode for an
        // integration (two copies per cell for each human chromosome pair),
        // the mean for DNA found and the phiX control; none otherwise.
        let copies = match r.call {
            ViralCall::Ici => Some(r.mode as f64),
            ViralCall::Found | ViralCall::LabControl => Some(r.mean),
            _ => None,
        };
        let mut cells = r.cells.clone();
        cells.push(
            copies
                .zip(hdepth)
                .map_or(".".into(), |(d, h)| crate::fixed(2.0 * d / h, 4)),
        );
        cells.push(r.call.text().into());
        cells.push(r.note.clone());
        w.write_row(cells).map_err(io)?;
    }
    crate::write_atomic(f.out, &w.into_inner())?;
    if !assessed {
        return Ok(format!(
            "{sample}: {} viral genomes; not assessed: the BAM has no unplaced unmapped read\n",
            rows.len()
        ));
    }
    let found: Vec<String> = rows
        .iter()
        .filter(|r| r.call != ViralCall::NotFound)
        .map(|r| match r.call {
            ViralCall::BelowThreshold => {
                format!("{} (below threshold, {} reads)", r.virus, r.reads)
            }
            _ => format!("{} ({})", r.virus, r.call.text()),
        })
        .collect();
    Ok(format!(
        "{sample}: {} viral genomes; {}\n",
        rows.len(),
        if found.is_empty() {
            "no reads above the threshold at this depth".to_owned()
        } else {
            found.join("; ")
        }
    ))
}

const SPEC: Spec = Spec {
    prog: "og viral report",
    usage: "usage: og viral report [-h] --sample SAMPLE --viral-fasta VIRAL_FASTA --viral-coverage VIRAL_COVERAGE --viral-depth VIRAL_DEPTH --human-coverage HUMAN_COVERAGE --viral-hits VIRAL_HITS --idxstats IDXSTATS --out OUT\n",
    options: &[
        "--sample",
        "--viral-fasta",
        "--viral-coverage",
        "--viral-depth",
        "--human-coverage",
        "--viral-hits",
        "--idxstats",
        "--out",
    ],
    positionals: &[],
    required: &[
        "--sample",
        "--viral-fasta",
        "--viral-coverage",
        "--viral-depth",
        "--human-coverage",
        "--viral-hits",
        "--idxstats",
        "--out",
    ],
};

/// `og viral VERB ARGS`: the exit code.
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if verb != Some("report") {
        return crate::unknown_verb("viral", "report", verb, err);
    }
    crate::run(&SPEC, rest, out, err, |a| {
        let files = Files {
            viral_fasta: crate::opt(a, "--viral-fasta"),
            viral_coverage: crate::opt(a, "--viral-coverage"),
            viral_depth: crate::opt(a, "--viral-depth"),
            human_coverage: crate::opt(a, "--human-coverage"),
            viral_hits: crate::opt(a, "--viral-hits"),
            idxstats: crate::opt(a, "--idxstats"),
            out: crate::opt(a, "--out"),
        };
        Ok(report(crate::opt(a, "--sample"), &files).map_err(|e| (e, 1)))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titles_and_names() {
        assert_eq!(
            short_title("Human gammaherpesvirus 4, complete genome"),
            "Human gammaherpesvirus 4"
        );
        assert_eq!(
            short_title("Hepatitis B virus genome sequence"),
            "Hepatitis B virus"
        );
        assert_eq!(
            short_title("Torque teno virus DNA, isolate x"),
            "Torque teno virus"
        );
        assert_eq!(short_title("Some genomes"), "Some genomes");
        assert_eq!(short_name("NC_001664.4", "x"), "HHV-6A");
        assert_eq!(short_name("NC_1", ", complete genome"), "NC_1");
        assert!(autosome("chr22") && autosome("7") && !autosome("chr23") && !autosome("chr0"));
    }

    #[test]
    fn mode_tie_goes_to_the_lower_depth() {
        let d = crate::testutil::Dir::new("viral");
        let p = d.put("d.tsv", "V\t1\t3\nV\t2\t3\nV\t3\t1\nV\t4\t1\n");
        let lengths: IndexMap<String, i128> = [("V".to_owned(), 4)].into_iter().collect();
        let got = read_depth(&p, &lengths).unwrap();
        assert_eq!(got["V"], (2.0, 1));
    }

    /// Evidence of `reads` counted reads over `covbases` of a `length`-base genome.
    fn ev(reads: i128, covbases: i128, length: i128, relative: Option<f64>) -> Evidence {
        Evidence {
            aligned: reads,
            counted: Counted { reads, covbases },
            length,
            relative,
        }
    }

    /// Red first: a relative depth of exactly 0.30 is integrated HHV-6.
    #[test]
    fn red_first_the_ici_threshold_is_inclusive() {
        assert_eq!(
            call_row("NC_000898.1", &ev(0, 0, 2000, Some(0.30))).0,
            ViralCall::Ici
        );
        assert_eq!(
            call_row("NC_000898.1", &ev(50, 1000, 2000, Some(0.2999))).0,
            ViralCall::Found
        );
        assert_eq!(
            call_row("NC_001422.1", &ev(3, 0, 500, None)).0,
            ViralCall::LabControl
        );
    }

    /// Red first (#411 trait-54): 1 to 9 counted reads read "below threshold"
    /// with their count, never "not found"; none at all is "not found".
    #[test]
    fn red_first_one_to_nine_reads_are_below_threshold() {
        let (call, note) = call_row("NC_007605.1", &ev(4, 600, 171_823, Some(0.0)));
        assert_eq!(call, ViralCall::BelowThreshold);
        assert_eq!(call.text(), "below threshold");
        assert!(note.starts_with("4 counted reads, under the 10"), "{note}");
        assert_eq!(
            call_row("NC_007605.1", &ev(9, 1500, 1500, Some(0.0))).0,
            ViralCall::BelowThreshold
        );
        assert_eq!(
            call_row("NC_007605.1", &ev(0, 0, 1500, Some(0.0))).0,
            ViralCall::NotFound
        );
    }

    /// Red first (#411 trait-52): exactly a tenth of the genome is "found";
    /// `1.0 - 135000.0 / 150000.0` is below 0.10 in floating point.
    #[test]
    fn red_first_a_tenth_of_the_genome_is_found() {
        assert_eq!(
            call_row("NC_006273.2", &ev(100, 15_000, 150_000, Some(0.0))).0,
            ViralCall::Found
        );
        assert_eq!(
            call_row("NC_006273.2", &ev(100, 14_999, 150_000, Some(0.0))).0,
            ViralCall::FewReads
        );
    }

    #[test]
    fn dust_scores_low_complexity_high() {
        let poly_a = vec![b'A'; 150];
        let ca: Vec<u8> = b"CA".iter().copied().cycle().take(150).collect();
        let telo: Vec<u8> = b"TTAGGG".iter().copied().cycle().take(150).collect();
        assert!((dust(&poly_a) - 100.0).abs() < 1e-9, "{}", dust(&poly_a));
        assert!(dust(&ca) > 40.0, "{}", dust(&ca));
        assert!(dust(&telo) > DUST_MAX, "{}", dust(&telo));
        // A fixed pseudo-random read stays far under the cut-off.
        let mut x = 47u32;
        let random: Vec<u8> = (0..150)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect();
        assert!(dust(&random) < 4.0, "{}", dust(&random));
        assert!(dust(b"").abs() < 1e-9 && dust(b"ACGTN").abs() < 1e-9);
    }

    #[test]
    fn sam_records_that_never_count() {
        let rec = |flag: u16, mapq: u8, cigar: &str, seq: &str| {
            format!("r\t{flag}\tV\t11\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*\n")
        };
        let seq = "ACGTTGCAAGCTTCGA";
        assert!(parse_hit(&rec(0, 60, "16M", seq)).unwrap().is_some());
        for (flag, mapq) in [(4, 60), (256, 60), (512, 60), (1024, 60), (0, 19), (0, 255)] {
            assert!(parse_hit(&rec(flag, mapq, "16M", seq)).unwrap().is_none());
        }
        assert!(
            parse_hit(&rec(0, 60, "16M", "AAAAAAAAAAAAAAAA"))
                .unwrap()
                .is_none()
        );
        // 5' ends: forward from POS minus the leading clip, reverse from the
        // end plus the trailing clip.
        let (_, rev, five, supp, blocks) =
            parse_hit(&rec(0, 60, "2S10M2D4M", seq)).unwrap().unwrap();
        assert_eq!(
            (rev, five, supp, blocks),
            (false, 9, false, vec![(11, 21), (23, 27)])
        );
        let (_, rev, five, _, _) = parse_hit(&rec(16, 60, "10M3H", seq)).unwrap().unwrap();
        assert_eq!((rev, five), (true, 23));
        assert!(parse_hit(&rec(2048, 60, "10M", seq)).unwrap().unwrap().3);
        assert!(parse_hit("r\t0\tV\n").is_err());
        assert!(parse_hit(&rec(0, 60, "10Q", seq)).is_err());
        assert!(parse_hit(&rec(0, 60, "5S", seq)).is_err());
    }
}
