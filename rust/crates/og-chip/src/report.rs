//! `og chip report`: each protein change of a table gene on its MANE
//! transcript, with its class and status, and per gene the depth and the
//! read-support limit.
//!
//! Inputs, written by bcftools (`scripts/49-chip.sh`, `CHIP_ANNOTATE`):
//! - calls: `CHROM POS REF ALT FILTER BCSQ AD AF DP`, one ALT per line,
//!   annotated by `bcftools csq --local-csq`;
//! - depth: `CHROM POS DP` from `bcftools mpileup` over the BED; a base with
//!   no line has depth 0.
//!
//! Statuses, in the order they are tested:
//! - `likely_germline`: VAF [`GERMLINE_VAF`] or more (an inherited variant
//!   sits near 50% or 100%, and a tumour-only call cannot tell it from a
//!   large clone);
//! - `filtered`: FilterMutectCalls did not pass it;
//! - `low_support`: fewer than [`MIN_ALT_READS`] reads with the change;
//! - `low_vaf`: VAF under [`CHIP_MIN_VAF`], the 2% that defines CHIP;
//! - `chip`: in the table's classes. The note says it is a possible finding,
//!   not a diagnosis, and adds a germline-fit sentence when its reads fit an
//!   inherited heterozygous variant ([`fits_inherited`]);
//! - `other_change`: a protein change outside them.
//!
//! A base counts once per gene, even when the depth file lists it twice.
//! Its depth goes into the mean of only the last gene of the BED that
//! covers it, as the Python's mean did, and into the per-base limit and
//! the shares of every gene that covers it, so a base two genes share is
//! read for both.
//!
//! Per gene, the limit comes from each coding base's own depth, not from
//! the mean ([`crate::binom::limit_over`]), so an exon the run did not read
//! raises the limit or leaves none:
//! - `covered`: a limit exists;
//! - `partly_covered`: some bases have [`MIN_ALT_READS`] reads, too few
//!   for any VAF to reach the power over the gene;
//! - `too_shallow`: no base has [`MIN_ALT_READS`] reads;
//! - `not_on_reference`: no coding base on the reference.
//!
//! A limit is written rounded up at the third decimal, so a clone at the
//! written VAF reaches the power. The stated limit of the run is the median
//! of the covered genes' written limits, the figure every reader of the
//! table takes, and the shares are counted at that same figure.
//! Each gene also gets the share of its coding bases with
//! [`CALLABLE_DEPTH`] reads or more, and the share where a clone at the
//! stated limit shows its reads.
//!
//! The calls table differs from the Python's on purpose since issue #411:
//! the 2% floor, the notes of `chip` and `likely_germline` rows, and the
//! ASXL1 row of the table. The depth table and the summary line differ on
//! purpose since the per-base limit.

use crate::binom::{detection_limit_over, finite_mean, readable_at};
use crate::prepare::{attr, columns};
use crate::table::{
    Missense, Row, Truncating, digits, parse_missense, parse_truncating, read_table,
};
use crate::{CALLABLE_DEPTH, DETECT_POWER, Error, MIN_ALT_READS, each_line, int, write_atomic};
use indexmap::IndexMap;
use pycompat::fmt::fixed;
use pycompat::num::{parse_float, py_strip};
use std::collections::{BTreeSet, HashMap};

/// At or above this VAF a call is read as a possible inherited variant.
pub const GERMLINE_VAF: f64 = 0.35;
/// CHIP is defined from a VAF of 2% (Steensma et al. 2015, PMID 25931582):
/// a call under it is `low_vaf`, not counted.
pub const CHIP_MIN_VAF: f64 = 0.02;
/// A `chip` call whose alt reads are this likely or more under an inherited
/// heterozygous variant (one-sided, `P(X <= alt)` with `X ~ Binomial(depth,
/// 0.5)`) gets the germline-fit note: 9 or 10 of 30 reads fit, 6 of 46 do
/// not.
pub const GERMLINE_FIT_P: f64 = 0.01;
/// What a `chip` row's note starts with: the step finds a possible clone,
/// it does not diagnose CHIP.
pub const CHIP_LABEL: &str = "possible clonal haematopoiesis (not a diagnosis)";

/// `P(X <= alt)` for `X ~ Binomial(depth, 0.5)`, summed in log space so a
/// deep base neither overflows nor underflows. 1 when `alt >= depth`.
pub fn het_lower_tail(alt: u64, depth: u64) -> f64 {
    if alt >= depth {
        return 1.0;
    }
    let n = depth as f64;
    // ln P(X = i), from ln P(X = 0) = n ln(1/2) by the ratio (n - i) / (i + 1).
    let mut l = n * 0.5f64.ln();
    let mut logs = Vec::with_capacity(alt as usize + 1);
    for i in 0..=alt {
        logs.push(l);
        l += ((n - i as f64) / (i as f64 + 1.0)).ln();
    }
    let top = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = logs.iter().map(|x| (x - top).exp()).sum();
    (top + sum.ln()).exp().min(1.0)
}

/// Whether `alt` of `depth` reads fit an inherited heterozygous variant:
/// [`het_lower_tail`] at [`GERMLINE_FIT_P`] or more.
pub fn fits_inherited(alt: u64, depth: u64) -> bool {
    depth > 0 && het_lower_tail(alt, depth) >= GERMLINE_FIT_P
}

/// The calls table's columns.
pub const OUT_COLUMNS: [&str; 15] = [
    "gene",
    "chrom",
    "pos",
    "ref",
    "alt",
    "consequence",
    "protein_change",
    "exon",
    "vaf",
    "alt_reads",
    "depth",
    "filter",
    "status",
    "class",
    "note",
];
/// The depth table's columns. The last two came after the first release
/// of the table, so a reader of the first seven still reads it.
pub const DEPTH_COLUMNS: [&str; 9] = [
    "gene",
    "transcript",
    "chrom",
    "coding_bases",
    "mean_depth",
    "detection_limit_vaf",
    "status",
    "callable_fraction",
    "assessable_fraction",
];

const TRUNCATING: [&str; 4] = [
    "frameshift",
    "stop_gained",
    "splice_donor",
    "splice_acceptor",
];
const PROTEIN: [&str; 11] = [
    "frameshift",
    "stop_gained",
    "splice_donor",
    "splice_acceptor",
    "missense",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_altering",
    "start_lost",
    "stop_lost",
    "feature_elongation",
];

/// A gene's MANE transcript, from prepare's GFF3.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GeneModel {
    pub tx: String,
    pub chrom: String,
    pub strand: String,
    /// `(start, end, rank)`, in file order.
    pub exons: Vec<(i128, i128, i128)>,
}

impl GeneModel {
    /// The rank of the exon holding `pos`, or of the nearest one (a splice
    /// site sits in the intron); 0 without exons.
    pub fn exon_of(&self, pos: i128) -> i128 {
        let mut best: Option<(i128, i128)> = None;
        for &(s, e, rank) in &self.exons {
            let d = if s <= pos && pos <= e {
                0
            } else {
                (pos - s).abs().min((pos - e).abs())
            };
            if best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, rank));
            }
        }
        best.map_or(0, |b| b.1)
    }
}

/// Gene symbol -> its transcript, in the GFF3's order.
pub fn read_gff(path: &str) -> Result<IndexMap<String, GeneModel>, Error> {
    let mut out: IndexMap<String, GeneModel> = IndexMap::new();
    let mut tx_gene: HashMap<String, String> = HashMap::new();
    each_line(path, |ln| {
        if ln.starts_with('#') {
            return Ok(());
        }
        let Some(f) = columns(ln.strip_suffix('\n').unwrap_or(ln)) else {
            return Ok(());
        };
        let c9 = f[8];
        let parent = attr(c9, "Parent").unwrap_or("");
        if f[2] == "gene" {
            tx_gene.insert(
                attr(c9, "ID").unwrap_or("").to_owned(),
                attr(c9, "Name").unwrap_or("").to_owned(),
            );
        } else if let Some(g) = tx_gene.get(parent)
            && let Some(tid) = attr(c9, "transcript_id")
        {
            let g = g.clone();
            out.insert(
                g.clone(),
                GeneModel {
                    tx: tid.to_owned(),
                    chrom: f[0].to_owned(),
                    strand: f[6].to_owned(),
                    exons: Vec::new(),
                },
            );
            tx_gene.insert(format!("transcript:{tid}"), g);
        } else if f[2] == "exon"
            && parent.starts_with("transcript:")
            && let Some(g) = tx_gene.get(parent)
            && !g.is_empty()
        {
            let rank = match attr(c9, "rank") {
                Some(r) => int(r, path)?,
                None => 0,
            };
            let exon = (int(f[3], path)?, int(f[4], path)?, rank);
            if let Some(m) = out.get_mut(g) {
                m.exons.push(exon);
            }
        }
        Ok(())
    })?;
    Ok(out)
}

/// What BCSQ says of one transcript.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Bcsq {
    pub csq: Vec<String>,
    /// `R882H`, `Q1200fs`, or the raw field when it does not parse.
    pub change: String,
    pub codon: i128,
    pub ref_aa: String,
    pub alt_aa: String,
}

/// `[A-Z*]*` from the start of `s`: its length.
fn residues(s: &[u8]) -> usize {
    s.iter()
        .take_while(|c| c.is_ascii_uppercase() || **c == b'*')
        .count()
}

fn digit_run(s: &[u8]) -> usize {
    s.iter().take_while(|c| c.is_ascii_digit()).count()
}

/// `(\d+)([A-Z*]*)(?:\.\.\d+[A-Z*]*)?>(\d+)?([A-Z*]*)` at the start of the
/// amino-acid field: codon digits, the reference residues and the new ones.
fn parse_aa(aa: &str) -> Option<(&str, &str, &str)> {
    let b = aa.as_bytes();
    let n = digit_run(b);
    if n == 0 {
        return None;
    }
    let r = residues(&b[n..]);
    let mut i = n + r;
    if b[i..].starts_with(b"..") {
        let d = digit_run(&b[i + 2..]);
        if d == 0 {
            return None;
        }
        i += 2 + d;
        i += residues(&b[i..]);
    }
    if b.get(i) != Some(&b'>') {
        return None;
    }
    i += 1;
    i += digit_run(&b[i..]);
    let a = residues(&b[i..]);
    Some((&aa[..n], &aa[n..n + r], &aa[i..i + a]))
}

/// The consequences and protein change of transcript `tx` in a BCSQ field
/// (`missense|DNMT3A|ENST00000321117|protein_coding|-|882R>882H|25234373C>T`
/// -> `missense`, `R882H`). A frameshift reads `R882fs`; a change of several
/// residues keeps the first one. A consequence marked downstream of a stop
/// (a leading `*`) or a compound pointer (`@`) is skipped.
pub fn parse_bcsq(field: &str, tx: &str) -> Bcsq {
    for entry in field.split(',') {
        let p: Vec<&str> = entry.split('|').collect();
        if p.len() < 3
            || p[2].split('.').next() != Some(tx)
            || entry.starts_with('@')
            || entry.starts_with('*')
        {
            continue;
        }
        let csq: Vec<String> = p[0].split('&').map(str::to_owned).collect();
        let aa = p.get(5).copied().unwrap_or("");
        return match parse_aa(aa) {
            Some((num, refs, alts)) => {
                let ref_aa = refs.get(..1).unwrap_or("").to_owned();
                let alt_aa = if refs.len() == 1 {
                    alts.to_owned()
                } else {
                    String::new()
                };
                let change = if csq.iter().any(|c| c == "frameshift") {
                    format!("{ref_aa}{num}fs")
                } else {
                    format!("{ref_aa}{num}{alt_aa}")
                };
                Bcsq {
                    csq,
                    change,
                    codon: digits(num).unwrap_or(0),
                    ref_aa,
                    alt_aa,
                }
            }
            None => Bcsq {
                csq,
                change: aa.to_owned(),
                ..Bcsq::default()
            },
        };
    }
    Bcsq::default()
}

/// The class a change belongs to (`""` when none) and why.
pub fn in_class(row: &Row, b: &Bcsq, exon: i128, nexons: usize) -> (&'static str, String) {
    let trunc = parse_truncating(&row.truncating).unwrap_or(Truncating::None);
    let mut hit: Vec<&str> = TRUNCATING
        .iter()
        .copied()
        .filter(|t| b.csq.iter().any(|c| c == t))
        .collect();
    hit.sort_unstable();
    if let Some(first) = hit.first() {
        return match trunc {
            Truncating::All => ("truncating", (*first).to_owned()),
            Truncating::Exons(a, z) if a <= exon && exon <= z => {
                ("truncating", format!("{first} in exon {exon}"))
            }
            Truncating::Exons(a, z) => (
                "",
                format!(
                    "truncating, but in exon {exon} of {nexons} (the table counts exons {a}-{z})"
                ),
            ),
            Truncating::None => (
                "",
                "truncating changes of this gene are not counted".to_owned(),
            ),
        };
    }
    if b.csq.iter().any(|c| c == "missense") {
        let one = |s: &str| s.as_bytes().first().copied();
        for rule in parse_missense(&row.missense).unwrap_or_default() {
            match rule {
                Missense::All => return ("missense", "missense".to_owned()),
                Missense::Range(a, z) if a <= b.codon && b.codon <= z => {
                    return ("missense", format!("missense in codons {a}-{z}"));
                }
                Missense::Exons(a, z) if a <= exon && exon <= z => {
                    return ("missense", format!("missense in exon {exon}"));
                }
                Missense::Codon(r, n)
                    if b.codon == n && b.ref_aa.len() == 1 && one(&b.ref_aa) == Some(r) =>
                {
                    return ("missense", format!("hotspot {}{n}", r as char));
                }
                Missense::Change(r, n, a)
                    if b.codon == n
                        && b.ref_aa.len() == 1
                        && one(&b.ref_aa) == Some(r)
                        && b.alt_aa.len() == 1
                        && one(&b.alt_aa) == Some(a) =>
                {
                    return ("missense", format!("{}{n}{}", r as char, a as char));
                }
                _ => {}
            }
        }
        return (
            "",
            "a missense change outside the table's hotspots".to_owned(),
        );
    }
    (
        "",
        "a protein change outside the table's classes".to_owned(),
    )
}

fn num(x: &str) -> Option<f64> {
    parse_float(x).ok()
}

/// `int(x)` of a float: truncated; a value that is not finite is an error.
fn float_int(x: f64, what: &str) -> Result<i128, Error> {
    if x.is_finite() {
        Ok(x.trunc() as i128)
    } else {
        Err(Error(format!(
            "{what} {} is not a number of reads",
            pycompat::fmt::repr(x)
        )))
    }
}

/// One row of the calls table, in [`OUT_COLUMNS`] order.
pub type OutRow = [String; 15];

/// The rows of the calls (`CHROM POS REF ALT FILTER BCSQ AD AF DP`), in
/// input order: one per call that changes a protein of a table gene.
pub fn classify(
    calls: &[Vec<String>],
    table: &[Row],
    genes: &IndexMap<String, GeneModel>,
    path: &str,
) -> Result<Vec<OutRow>, Error> {
    let rows: HashMap<&str, &Row> = table.iter().map(|r| (r.gene.as_str(), r)).collect();
    let mut out = Vec::new();
    for c in calls {
        let field = |i: usize| c.get(i).map_or("", String::as_str);
        let (chrom, pos, flt, ad, af, dp) =
            (field(0), field(1), field(4), field(6), field(7), field(8));
        for (g, info) in genes {
            let Some(row) = rows.get(g.as_str()) else {
                continue;
            };
            if info.chrom != chrom {
                continue;
            }
            let b = parse_bcsq(field(5), &info.tx);
            if !b.csq.iter().any(|c| PROTEIN.contains(&c.as_str())) {
                continue;
            }
            let exon = info.exon_of(int(pos, path)?);
            let (cls, why) = in_class(row, &b, exon, info.exons.len());
            let ads: Vec<Option<f64>> = if ad.is_empty() || ad == "." {
                Vec::new()
            } else {
                ad.split(',').map(num).collect()
            };
            let alt_reads = match ads.get(1) {
                Some(Some(x)) => float_int(*x, "AD")?,
                _ => 0,
            };
            let mut depth = num(dp);
            if depth.is_none() && !ads.is_empty() && ads.iter().all(Option::is_some) {
                depth = Some(pycompat::num::sum(ads.iter().flatten().copied()));
            }
            let mut vaf = if af.is_empty() || af == "." {
                None
            } else {
                num(af.split(',').next().unwrap_or(""))
            };
            if vaf.is_none()
                && let Some(d) = depth
                && d != 0.0
            {
                vaf = Some(alt_reads as f64 / d);
            }
            let (status, note) = if let Some(v) = vaf.filter(|v| *v >= GERMLINE_VAF) {
                (
                    "likely_germline",
                    format!(
                        "VAF {}: at or above the 35% cut, where an inherited variant (near 50% or 100%) and a large clone look the same without a matched normal sample; not counted",
                        fixed(v, 2)
                    ),
                )
            } else if flt != "PASS" && flt != "." {
                ("filtered", format!("FilterMutectCalls: {flt}"))
            } else if alt_reads < i128::from(MIN_ALT_READS) {
                (
                    "low_support",
                    format!("{alt_reads} reads with the change, under {MIN_ALT_READS}"),
                )
            } else if let Some(v) = vaf.filter(|v| *v < CHIP_MIN_VAF) {
                (
                    "low_vaf",
                    format!(
                        "VAF {}: under the 2% VAF that defines CHIP, not counted",
                        fixed(v, 3)
                    ),
                )
            } else if !cls.is_empty() {
                let mut note = format!("{CHIP_LABEL}: {why}");
                // The reads the table shows: alt_reads of the depth column.
                if let Some(d) = depth.filter(|d| d.is_finite() && *d >= 1.0)
                    && alt_reads >= 0
                    && fits_inherited(alt_reads as u64, d as u64)
                {
                    note.push_str(&format!(
                        "; {alt_reads} of {} reads also fit an inherited heterozygous variant",
                        d as u64
                    ));
                }
                ("chip", note)
            } else {
                ("other_change", why)
            };
            let depth = match depth {
                Some(d) => float_int(d, "DP")?.to_string(),
                None => String::new(),
            };
            out.push([
                g.clone(),
                chrom.to_owned(),
                pos.to_owned(),
                field(2).to_owned(),
                field(3).to_owned(),
                b.csq.join("&"),
                b.change,
                exon.to_string(),
                vaf.map_or(String::new(), |v| fixed(v, 3)),
                alt_reads.to_string(),
                depth,
                flt.to_owned(),
                status.to_owned(),
                cls.to_owned(),
                note,
            ]);
            break;
        }
    }
    Ok(out)
}

/// The bases of one contig the BED covers, cut where the set of covering
/// rows changes: `(first, last + 1, owner, genes)`, 1-based, ascending and
/// disjoint, where `owner` is the gene whose mean counts the base and
/// `genes` every gene that covers it, with a flag per base for "its depth
/// was counted".
struct Segments {
    segs: Vec<(i128, i128, usize, Vec<usize>)>,
    counted: Vec<Vec<bool>>,
}

impl Segments {
    /// `rows`: `(first, last + 1, gene index, order)`; where rows overlap,
    /// the one with the highest order owns the base.
    fn new(mut rows: Vec<(i128, i128, usize, usize)>) -> Segments {
        rows.retain(|r| r.0 < r.1);
        let mut cuts: Vec<i128> = rows.iter().flat_map(|r| [r.0, r.1]).collect();
        cuts.sort_unstable();
        cuts.dedup();
        let mut starts = rows.clone();
        starts.sort_by_key(|r| r.0);
        let mut ends = rows;
        ends.sort_by_key(|r| r.1);
        let (mut si, mut ei) = (0, 0);
        let mut active: BTreeSet<(usize, usize)> = BTreeSet::new();
        let mut segs = Vec::new();
        for w in cuts.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            while ei < ends.len() && ends[ei].1 <= lo {
                active.remove(&(ends[ei].3, ends[ei].2));
                ei += 1;
            }
            while si < starts.len() && starts[si].0 <= lo {
                active.insert((starts[si].3, starts[si].2));
                si += 1;
            }
            if let Some(&(_, owner)) = active.last() {
                let mut genes: Vec<usize> = active.iter().map(|&(_, g)| g).collect();
                genes.sort_unstable();
                genes.dedup();
                segs.push((lo, hi, owner, genes));
            }
        }
        let counted = segs
            .iter()
            .map(|s| vec![false; (s.1 - s.0) as usize])
            .collect();
        Segments { segs, counted }
    }

    /// The owner of `pos`, every gene that covers it, and its counted flag.
    fn at(&mut self, pos: i128) -> Option<(usize, &[usize], &mut bool)> {
        let i = self.segs.partition_point(|s| s.1 <= pos);
        let (lo, hi, owner, genes) = self.segs.get(i)?;
        if pos < *lo || pos >= *hi {
            return None;
        }
        Some((*owner, genes, &mut self.counted[i][(pos - lo) as usize]))
    }
}

/// Bases covered by half-open `(start, end)` rows.
fn bases(mut rows: Vec<(i128, i128)>) -> i128 {
    rows.retain(|r| r.0 < r.1);
    rows.sort_unstable();
    let (mut total, mut cur): (i128, Option<(i128, i128)>) = (0, None);
    for (s, e) in rows {
        match cur {
            Some((cs, ce)) if s <= ce => cur = Some((cs, ce.max(e))),
            Some((cs, ce)) => {
                total += ce - cs;
                cur = Some((s, e));
            }
            None => cur = Some((s, e)),
        }
    }
    total + cur.map_or(0, |(s, e)| e - s)
}

/// One row of the depth table, in [`DEPTH_COLUMNS`] order.
pub type DepthRow = [String; 9];

/// What one gene's coding bases hold, for its depth row.
struct GeneDepth {
    /// Coding bases on the reference.
    n: u64,
    mean: f64,
    /// `(depth, bases)`, the unread bases at depth 0.
    depths: Vec<(u64, u64)>,
    limit: Option<f64>,
}

impl GeneDepth {
    /// The share of the coding bases where `keep` holds for the depth.
    fn share(&self, keep: impl Fn(u64) -> bool) -> String {
        let k: u64 = self
            .depths
            .iter()
            .filter(|&&(d, _)| keep(d))
            .map(|&(_, b)| b)
            .sum();
        fixed(k as f64 / self.n as f64, 3)
    }
}

/// A depth as whole reads: a fraction is cut, a negative or `nan` is 0.
fn reads(d: f64) -> u64 {
    // `as` truncates, saturates and turns NaN into 0.
    d as u64
}

/// One row per table gene: its coding bases on the reference, the mean
/// depth over them, the read-support limit from each base's depth, its
/// status, and the shares of its bases at [`CALLABLE_DEPTH`] reads and at
/// the stated limit of the run.
pub fn depth_rows(
    table: &[Row],
    genes: &IndexMap<String, GeneModel>,
    bed: &str,
    depth: &str,
) -> Result<Vec<DepthRow>, Error> {
    // Gene -> its BED rows (contig, start, end), genes in first-seen order.
    let mut regions: IndexMap<String, Vec<(String, i128, i128)>> = IndexMap::new();
    each_line(bed, |ln| {
        let f: Vec<&str> = ln.strip_suffix('\n').unwrap_or(ln).split('\t').collect();
        if f.len() >= 4 {
            let row = (f[0].to_owned(), int(f[1], bed)?, int(f[2], bed)?);
            regions.entry(f[3].to_owned()).or_default().push(row);
        }
        Ok(())
    })?;
    // Per contig, who covers each base, and who owns it for the mean: the
    // last assignment in gene order, then row order, as the Python's
    // position dictionary kept it.
    let mut by_chrom: HashMap<&str, Vec<(i128, i128, usize, usize)>> = HashMap::new();
    let mut order = 0;
    for (gi, rows) in regions.values().enumerate() {
        for (c, s, e) in rows {
            by_chrom
                .entry(c)
                .or_default()
                .push((s + 1, e + 1, gi, order));
            order += 1;
        }
    }
    let mut owner: HashMap<String, Segments> = by_chrom
        .into_iter()
        .map(|(c, rows)| (c.to_owned(), Segments::new(rows)))
        .collect();
    let names: Vec<&str> = regions.keys().map(String::as_str).collect();
    let mut total = vec![0.0f64; regions.len()];
    let mut hist: Vec<HashMap<u64, u64>> = vec![HashMap::new(); regions.len()];
    each_line(depth, |ln| {
        let f: Vec<&str> = ln.strip_suffix('\n').unwrap_or(ln).split('\t').collect();
        if f.len() < 3 {
            return Ok(());
        }
        let pos = int(f[1], depth)?;
        if let Some(seg) = owner.get_mut(f[0])
            && let Some((gi, genes, counted)) = seg.at(pos)
            && !*counted
            && let Some(d) = num(f[2])
        {
            if !names[gi].is_empty() {
                total[gi] += d;
            }
            for &g in genes.iter().filter(|&&g| !names[g].is_empty()) {
                *hist[g].entry(reads(d)).or_default() += 1;
            }
            *counted = true;
        }
        Ok(())
    })?;

    let mut per_gene: Vec<(&Row, Option<GeneDepth>)> = Vec::with_capacity(table.len());
    for r in table {
        let gi = regions.get_index_of(r.gene.as_str());
        let n = gi.map_or(0, |i| {
            let rows = &regions[i];
            let mut per: HashMap<&str, Vec<(i128, i128)>> = HashMap::new();
            for (c, s, e) in rows {
                per.entry(c).or_default().push((*s, *e));
            }
            per.into_values().map(bases).sum()
        });
        if n <= 0 {
            per_gene.push((r, None));
            continue;
        }
        let n = n as u64;
        let mean = gi.map_or(0.0, |i| total[i]) / n as f64;
        finite_mean(mean)?;
        let mut depths: Vec<(u64, u64)> = gi
            .map(|i| hist[i].iter().map(|(&d, &b)| (d, b)).collect())
            .unwrap_or_default();
        let read: u64 = depths.iter().map(|&(_, b)| b).sum();
        if read < n {
            depths.push((0, n - read));
        }
        depths.sort_unstable();
        let limit = detection_limit_over(&depths);
        per_gene.push((
            r,
            Some(GeneDepth {
                n,
                mean,
                depths,
                limit,
            }),
        ));
    }
    let stated = stated_limit(
        per_gene
            .iter()
            .filter_map(|(_, g)| g.as_ref()?.limit.map(written_limit)),
    );

    let mut out = Vec::with_capacity(per_gene.len());
    for (r, gd) in per_gene {
        let g = r.gene.as_str();
        let info = genes.get(g).cloned().unwrap_or_default();
        let Some(gd) = gd else {
            out.push([
                g.to_owned(),
                info.tx,
                info.chrom,
                "0".to_owned(),
                String::new(),
                String::new(),
                "not_on_reference".to_owned(),
                String::new(),
                String::new(),
            ]);
            continue;
        };
        let status = if gd.limit.is_some() {
            "covered"
        } else if gd
            .depths
            .iter()
            .any(|&(d, b)| b > 0 && d >= u64::from(MIN_ALT_READS))
        {
            "partly_covered"
        } else {
            "too_shallow"
        };
        out.push([
            g.to_owned(),
            info.tx,
            info.chrom,
            gd.n.to_string(),
            fixed(gd.mean, 1),
            gd.limit
                .map_or(String::new(), |l| fixed(written_limit(l), 3)),
            status.to_owned(),
            gd.share(|d| d >= CALLABLE_DEPTH),
            stated.map_or(String::new(), |f| gd.share(|d| readable_at(d, f))),
        ]);
    }
    Ok(out)
}

/// A limit as the depth table writes it: rounded up at the third decimal,
/// so the power at the written VAF is [`DETECT_POWER`] or more. The `1e-9`
/// keeps a limit that is already a whole thousandth from moving up by one.
fn written_limit(l: f64) -> f64 {
    (l * 1000.0 - 1e-9).ceil() / 1000.0
}

/// The stated limit of a run: the median of the covered genes' limits
/// (the mean of the middle two for an even count); `None` with none.
pub fn stated_limit(limits: impl IntoIterator<Item = f64>) -> Option<f64> {
    pycompat::stats::median(limits).ok()
}

/// The read-support sentence for the depth rows, read back from the rows
/// as `og summary collect` reads them.
pub fn summary_line(depth: &[DepthRow]) -> String {
    let covered: Vec<&DepthRow> = depth.iter().filter(|r| r[6] == "covered").collect();
    let (Some(med), Some(lim)) = (
        stated_limit(covered.iter().filter_map(|r| num(&r[4]))),
        stated_limit(covered.iter().filter_map(|r| num(&r[5]))),
    ) else {
        return format!(
            "No CHIP gene has the depth to find a clone (at least {MIN_ALT_READS} reads with the change): no result either way."
        );
    };
    let short: Vec<String> = depth
        .iter()
        .filter(|r| r[6] == "partly_covered" || r[6] == "too_shallow")
        .map(|r| match num(&r[8]) {
            Some(a) => format!("{} {}%", r[0], fixed(a * 100.0, 0)),
            None => r[0].clone(),
        })
        .collect();
    let short = if short.is_empty() {
        String::new()
    } else {
        format!(
            " Too few reads for any limit in {} gene{}; the share of each one's coding bases read deep enough for that VAF: {}.",
            short.len(),
            if short.len() == 1 { "" } else { "s" },
            short.join(", ")
        )
    };
    format!(
        "Read-support limit (median {}x over the CHIP genes' coding exons): in the median CHIP gene, a clone with a VAF of about {}% or more shows at least {MIN_ALT_READS} reads with the change {} times in 10, counted base by base. Mutect2 and its filters call fewer than that.{short} CHIP starts at a VAF of {}%, so smaller clones are likely missed: no call here does not rule one out. Not a diagnosis.",
        fixed(med, 0),
        fixed(lim * 100.0, 0),
        (DETECT_POWER * 10.0) as i64,
        fixed(CHIP_MIN_VAF * 100.0, 0)
    )
}

fn tsv<const N: usize>(columns: &[&str; N], rows: &[[String; N]]) -> String {
    let mut s = columns.join("\t");
    s.push('\n');
    for r in rows {
        let cells: Vec<String> = r.iter().map(|c| c.replace('\t', " ")).collect();
        s.push_str(&cells.join("\t"));
        s.push('\n');
    }
    s
}

/// The files `og chip report` reads and writes.
pub struct Files {
    pub genes: String,
    pub gff: String,
    pub bed: String,
    pub calls: String,
    pub depth: String,
    pub out: String,
    pub depth_out: String,
}

/// Lines with text, split on tabs, as the Python's `read_lines` read them.
fn read_calls(path: &str) -> Result<Vec<Vec<String>>, Error> {
    let mut out = Vec::new();
    each_line(path, |ln| {
        if !py_strip(ln).is_empty() {
            out.push(
                ln.strip_suffix('\n')
                    .unwrap_or(ln)
                    .split('\t')
                    .map(str::to_owned)
                    .collect(),
            );
        }
        Ok(())
    })?;
    Ok(out)
}

/// `og chip report`: writes both tables (each whole or not at all) and
/// returns the two lines for stdout.
pub fn report(files: &Files, sample: &str) -> Result<String, Error> {
    let table = read_table(&files.genes)?;
    let genes = read_gff(&files.gff)?;
    let rows = classify(&read_calls(&files.calls)?, &table, &genes, &files.calls)?;
    let depth = depth_rows(&table, &genes, &files.bed, &files.depth)?;
    let summary = summary_line(&depth);
    write_atomic(&files.out, tsv(&OUT_COLUMNS, &rows).as_bytes())?;
    write_atomic(&files.depth_out, tsv(&DEPTH_COLUMNS, &depth).as_bytes())?;
    let n = rows.iter().filter(|r| r[12] == "chip").count();
    let fit = rows
        .iter()
        .filter(|r| r[12] == "chip" && r[14].contains("also fit an inherited"))
        .count();
    let other = rows.len() - n;
    let s = |k: usize| if k == 1 { "" } else { "s" };
    let fit = if fit == 0 {
        String::new()
    } else {
        format!(
            "{fit} of them also fit{} an inherited variant; ",
            if fit == 1 { "s" } else { "" }
        )
    };
    Ok(format!(
        "{sample}: {n} CHIP call{} in {} genes (possible clonal haematopoiesis, not a diagnosis; {fit}{other} other protein change{} listed, not counted)\n{summary}\n",
        s(n),
        table.len(),
        s(other)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bcsq_of_one_transcript() {
        let b = parse_bcsq(
            "missense|DNMT3A|ENST00000321117|protein_coding|-|882R>882H|25234373C>T",
            "ENST00000321117",
        );
        assert_eq!(
            (
                b.change.as_str(),
                b.codon,
                b.ref_aa.as_str(),
                b.alt_aa.as_str()
            ),
            ("R882H", 882, "R", "H")
        );
        let fs = parse_bcsq("frameshift|TET2|T.3|protein_coding|+|1200QRS>1200P*|x", "T");
        assert_eq!((fs.change.as_str(), fs.alt_aa.as_str()), ("Q1200fs", ""));
        let span = parse_bcsq("inframe_deletion|G|T|pc|+|10AB..12C>10A|x", "T");
        assert_eq!(span.change, "A10");
        assert_eq!(
            parse_bcsq("*missense|G|T|pc|+|1R>1H|x", "T"),
            Bcsq::default()
        );
        assert_eq!(
            parse_bcsq("missense|G|U|pc|+|1R>1H|x", "T"),
            Bcsq::default()
        );
        let raw = parse_bcsq("splice_acceptor|G|T|pc|-", "T");
        assert_eq!(
            (raw.csq, raw.change, raw.codon),
            (vec!["splice_acceptor".to_owned()], String::new(), 0)
        );
        assert_eq!(
            parse_bcsq("missense|G|T|pc|+|1..x>1H|x", "T").change,
            "1..x>1H"
        );
    }

    #[test]
    fn nearest_exon() {
        let g = GeneModel {
            exons: vec![(100, 200, 2), (300, 400, 1)],
            ..GeneModel::default()
        };
        assert_eq!(g.exon_of(150), 2);
        assert_eq!(g.exon_of(250), 2);
        assert_eq!(g.exon_of(251), 1);
        assert_eq!(GeneModel::default().exon_of(5), 0);
    }

    #[test]
    fn a_base_belongs_to_the_last_row_over_it() {
        let mut s = Segments::new(vec![(1, 11, 0, 0), (6, 21, 1, 1), (3, 5, 2, 2)]);
        assert_eq!(s.at(1).map(|x| x.0), Some(0));
        assert_eq!(s.at(4).map(|x| x.0), Some(2));
        assert_eq!(s.at(7).map(|x| x.0), Some(1));
        assert_eq!(s.at(20).map(|x| x.0), Some(1));
        assert!(s.at(21).is_none() && s.at(0).is_none());
        assert_eq!(bases(vec![(0, 10), (5, 20), (30, 31), (40, 40)]), 21);
    }

    #[test]
    fn red_first_a_vaf_of_0_35_reads_as_germline() {
        let mut genes = IndexMap::new();
        genes.insert(
            "G".to_owned(),
            GeneModel {
                tx: "T".into(),
                chrom: "chr2".into(),
                strand: "-".into(),
                exons: vec![(100, 200, 1)],
            },
        );
        let at = |af: &str| {
            let c: Vec<String> = [
                "chr2",
                "150",
                "C",
                "T",
                "PASS",
                "missense|G|T|pc|-|882R>882H|x",
                "30,15",
                af,
                "45",
            ]
            .map(str::to_owned)
            .to_vec();
            classify(&[c], &[row("none", "R882")], &genes, "calls.tsv").unwrap()[0][12].clone()
        };
        assert_eq!(at("0.34"), "chip");
        assert_eq!(at("0.35"), "likely_germline");
        assert_eq!(at("0.97"), "likely_germline");
    }

    fn row(trunc: &str, mis: &str) -> Row {
        Row {
            gene: "G".into(),
            truncating: trunc.into(),
            missense: mis.into(),
            note: "n".into(),
            pmids: "123456".into(),
        }
    }

    #[test]
    fn classes() {
        let b = |csq: &str, aa: &str| parse_bcsq(&format!("{csq}|G|T|pc|+|{aa}|x"), "T");
        assert_eq!(
            in_class(
                &row("all", "none"),
                &b("stop_gained&splice_donor", "1R>1*"),
                3,
                5
            )
            .1,
            "splice_donor"
        );
        assert_eq!(
            in_class(&row("exons:6-6", "none"), &b("stop_gained", "1R>1*"), 5, 6).0,
            ""
        );
        assert_eq!(
            in_class(&row("none", "R882"), &b("missense", "882R>882H"), 1, 1).1,
            "hotspot R882"
        );
        assert_eq!(
            in_class(&row("none", "V617F"), &b("missense", "617V>617L"), 1, 1).0,
            ""
        );
        assert_eq!(
            in_class(&row("none", "102-292"), &b("missense", "248R>248Q"), 1, 1).1,
            "missense in codons 102-292"
        );
    }
}
