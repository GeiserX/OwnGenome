//! The depth-ratio engine: the guards in their order, the estimate, its 95%
//! interval and the call (`copy_number.assess`), and the PSV sites of a
//! `call psv` row (`copy_number.psv_assess`). The arithmetic is the
//! Python's step for step, with CPython 3.12's compensated `sum()` for
//! every interpreter (Parascopy's Python 3.11 summed plainly; the two agree
//! at every printed precision).

use crate::catalog::{CallRule, ControlMode, Locus, Segment, TargetMode};
use crate::round_i;
use pycompat::num::sum as py_sum;
use std::collections::HashMap;
use std::fmt;

/// Mean depth per segment name.
pub type Means = HashMap<String, f64>;

/// z of a two-sided 95% interval.
pub const Z95: f64 = 1.96;
/// Reads a PSV site needs (gene and pseudogene bases pooled) before its
/// copies are read.
pub const PSV_MIN_READS: u64 = 30;
/// A site reads as converted when its copies estimate lies this far under
/// the others' count (or under the usual count) or more. At 30x, about 70
/// reads at a site, sampling alone gives an estimate of 2 copies a standard
/// deviation of about 0.25; HG002's six sites read 1.62 to 2.82 (35x), so
/// 0.6 is a screen setting, not a measured error rate.
pub const PSV_MARGIN: f64 = 0.6;

/// The sample's sex, which sets the copies of chrX and chrY.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Sex {
    Male,
    Female,
    #[default]
    Unknown,
}

impl Sex {
    /// `male`, `female` or `unknown`.
    pub fn parse(s: &str) -> Option<Sex> {
        match s {
            "male" => Some(Sex::Male),
            "female" => Some(Sex::Female),
            "unknown" => Some(Sex::Unknown),
            _ => None,
        }
    }
}

/// What a block says.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Called,
    NotAssessed,
    /// A chrY row on a female (reason `no_y`).
    NotApplicable,
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Status::Called => "called",
            Status::NotAssessed => "not assessed",
            Status::NotApplicable => "not applicable",
        })
    }
}

/// Why a block is not called: the guard that stopped it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reason {
    NoData,
    LowControlDepth,
    ControlsDisagree,
    LowMapq0Share,
    HighMapq0Share,
    BelowFloor,
    BetweenCounts,
    SexUnknown,
    NoY,
    BrokenControl,
    CheckFailed,
    PsvLowDepth,
}

impl fmt::Display for Reason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Reason::NoData => "no_data",
            Reason::LowControlDepth => "low_control_depth",
            Reason::ControlsDisagree => "controls_disagree",
            Reason::LowMapq0Share => "low_mapq0_share",
            Reason::HighMapq0Share => "high_mapq0_share",
            Reason::BelowFloor => "below_floor",
            Reason::BetweenCounts => "between_counts",
            Reason::SexUnknown => "sex_unknown",
            Reason::NoY => "no_y",
            Reason::BrokenControl => "broken_control",
            Reason::CheckFailed => "check_failed",
            Reason::PsvLowDepth => "psv_low_depth",
        })
    }
}

/// A `call class` result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    Low,
    Typical,
    High,
}

impl Class {
    pub fn as_str(self) -> &'static str {
        match self {
            Class::Low => "low",
            Class::Typical => "typical",
            Class::High => "high",
        }
    }
}

/// What the PSV sites of a `call psv` row show.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Finding {
    None,
    Deletion,
    Conversion,
    PseudogeneDeletion,
}

impl Finding {
    pub fn as_str(self) -> &'static str {
        match self {
            Finding::None => "none",
            Finding::Deletion => "deletion",
            Finding::Conversion => "conversion",
            Finding::PseudogeneDeletion => "pseudogene_deletion",
        }
    }
}

/// One PSV site as read.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteRead {
    pub name: String,
    pub gene_reads: u64,
    pub pseudo_reads: u64,
    pub other_reads: u64,
    pub copies: Option<i64>,
    pub estimate: Option<f64>,
    pub ratio: Option<f64>,
}

/// A locus's result.
#[derive(Debug, Clone, PartialEq)]
pub struct Block {
    pub status: Status,
    pub reason: Option<Reason>,
    pub estimate: Option<f64>,
    pub total: Option<f64>,
    pub copies: Option<i64>,
    pub ci95: Option<(f64, f64)>,
    pub class: Option<Class>,
    pub target_depth: Option<f64>,
    pub control_depth: Option<f64>,
    /// Mean depth over every segment read (every read).
    pub depths: Means,
    /// The share of the target depth from reads of MAPQ 0; `None` without a
    /// MAPQ >= 1 pass, or when a guard stopped before the depths.
    pub mapq0_share: Option<f64>,
    pub reading: String,
    /// The copies with no change, when the sample's count of the chromosome
    /// is known.
    pub expected: Option<i64>,
    /// A `call psv` row: the gene and pseudogene copies together.
    pub module_total: Option<i64>,
    /// A `call psv` row: its sites, in catalog order (after `psv_assess`).
    pub sites: Vec<SiteRead>,
    pub finding: Option<Finding>,
    pub conversion_sites: Vec<String>,
}

impl Block {
    /// A block that says nothing yet (not assessed, no reason).
    pub fn empty() -> Block {
        Block {
            status: Status::NotAssessed,
            reason: None,
            estimate: None,
            total: None,
            copies: None,
            ci95: None,
            class: None,
            target_depth: None,
            control_depth: None,
            depths: Means::new(),
            mapq0_share: None,
            reading: String::new(),
            expected: None,
            module_total: None,
            sites: Vec::new(),
            finding: None,
            conversion_sites: Vec::new(),
        }
    }

    fn stop(mut self, reason: Reason) -> Block {
        self.reason = Some(reason);
        self
    }
}

/// Why `assess` gave no block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssessError {
    /// A control segment has no read while the other controls are covered:
    /// the row's coordinates are wrong, or the sample lacks that stretch.
    BrokenControl(String),
    /// A row with a MAPQ 0 guard was given no MAPQ >= 1 depths.
    NeedsMapq1(String),
}

impl fmt::Display for AssessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AssessError::BrokenControl(s) | AssessError::NeedsMapq1(s) => f.write_str(s),
        }
    }
}

fn combine_targets(mode: TargetMode, segs: &[Segment], m: &dyn Fn(&str) -> f64) -> f64 {
    match mode {
        TargetMode::Sum => py_sum(segs.iter().map(|s| m(&s.name))),
        TargetMode::Pooled => pooled(segs, m),
    }
}

fn combine_controls(mode: ControlMode, segs: &[&Segment], m: &dyn Fn(&str) -> f64) -> f64 {
    match mode {
        ControlMode::Flanks => py_sum(segs.iter().map(|s| m(&s.name))) / segs.len() as f64,
        ControlMode::Pooled => {
            let owned: Vec<Segment> = segs.iter().map(|s| (*s).clone()).collect();
            pooled(&owned, m)
        }
    }
}

/// Length-weighted mean: `sum(m * len) / sum(len)`.
fn pooled(segs: &[Segment], m: &dyn Fn(&str) -> f64) -> f64 {
    let num = py_sum(segs.iter().map(|s| m(&s.name) * (s.end - s.start) as f64));
    let den: i64 = segs.iter().map(|s| s.end - s.start).sum();
    num / den as f64
}

/// Coefficient of variation: sample standard deviation over the mean; 0
/// with fewer than two values or a mean of 0 or less.
fn cv(values: &[f64]) -> f64 {
    let mu = py_sum(values.iter().copied()) / values.len() as f64;
    if mu <= 0.0 || values.len() < 2 {
        return 0.0;
    }
    (py_sum(values.iter().map(|v| (v - mu).powi(2))) / (values.len() - 1) as f64).sqrt() / mu
}

/// The copies of the row's chromosome a sample has: 2, or 1 on a male's
/// chrX or chrY; `None` where it cannot be said (a chrY row on a female or
/// on a sample of no known sex, a chrX row without a sex).
pub fn ploidy(row: &Locus, sex: Sex) -> Option<i64> {
    let x = row.chrom == "chrX";
    if !x && row.chrom != "chrY" {
        return Some(2);
    }
    match sex {
        Sex::Male => Some(1),
        Sex::Female if x => Some(2),
        _ => None,
    }
}

/// The locus's block for mean depths `m` (every read) and `q1` (reads of
/// MAPQ >= 1, or `None` for a row without a MAPQ 0 guard).
pub fn assess(row: &Locus, m: &Means, q1: Option<&Means>, sex: Sex) -> Result<Block, AssessError> {
    let dp = row.depth();
    if q1.is_none() && (dp.mapq0_min.is_some() || dp.mapq0_max.is_some()) {
        return Err(AssessError::NeedsMapq1(format!(
            "{} has a MAPQ 0 guard: the MAPQ >= 1 depths are needed",
            row.id
        )));
    }
    let mut b = Block::empty();
    b.depths = m.clone();
    let units = if row.target_depth == "sum" {
        row.targets.len() as i64
    } else {
        1
    };
    let pl = ploidy(row, sex);
    b.expected = pl.map(|p| p * units);
    if row.chrom == "chrY" && sex == Sex::Female {
        b.status = Status::NotApplicable;
        return Ok(b.stop(Reason::NoY));
    }
    let Some(pl) = pl else {
        return Ok(b.stop(Reason::SexUnknown));
    };
    let (tg, ct) = (&row.targets, &row.controls);
    if m.is_empty() || tg.iter().chain(ct).any(|s| !m.contains_key(&s.name)) {
        return Ok(b.stop(Reason::NoData));
    }
    let get = |n: &str| m[n];
    for s in ct {
        let rest: Vec<&Segment> = ct.iter().filter(|c| c.name != s.name).collect();
        if m[&s.name] == 0.0 && !rest.is_empty() {
            let other = combine_controls(dp.control_mode, &rest, &get);
            if other >= dp.min_control && dp.min_control > 0.0 {
                return Err(AssessError::BrokenControl(format!(
                    "BROKEN_CONTROL: control segment {} of {} has no read while the other controls have {}x; check its coordinates in bin/cn_loci.tsv (or the sample lacks that stretch)",
                    s.name,
                    row.id,
                    crate::fixed(other, 1)
                )));
            }
        }
    }
    let all_ct: Vec<&Segment> = ct.iter().collect();
    let t = combine_targets(dp.target_mode, tg, &get);
    let c = combine_controls(dp.control_mode, &all_ct, &get);
    let share = q1.map(|q| {
        let tq = combine_targets(dp.target_mode, tg, &|n: &str| {
            q.get(n).copied().unwrap_or(0.0)
        });
        if t > 0.0 { 1.0 - tq / t } else { 0.0 }
    });
    b.target_depth = Some(t);
    b.control_depth = Some(c);
    b.mapq0_share = share;
    let low = match dp.control_mode {
        ControlMode::Flanks => ct.iter().map(|s| m[&s.name]).fold(f64::INFINITY, f64::min),
        ControlMode::Pooled => c,
    };
    if low < dp.min_control || c <= 0.0 {
        return Ok(b.stop(Reason::LowControlDepth));
    }
    if dp.control_mode == ControlMode::Flanks
        && ct.len() > 1
        && let Some(max_ratio) = dp.max_ratio
    {
        let vals: Vec<f64> = ct.iter().map(|s| m[&s.name]).collect();
        let lo = vals.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if lo <= 0.0 || hi / lo > max_ratio {
            return Ok(b.stop(Reason::ControlsDisagree));
        }
    }
    let total = dp.scale * pl as f64 / 2.0 * t / c;
    let est = total - dp.offset;
    b.total = Some(total);
    b.estimate = Some(est);
    // With no copy left there is no read to judge the aligner by (a row with
    // an offset has a floor of its own below, after this guard).
    let no_copy = dp.offset == 0.0 && round_i(est) == 0;
    if let (Some(min), Some(s)) = (dp.mapq0_min, share)
        && !no_copy
        && s < min
    {
        return Ok(b.stop(Reason::LowMapq0Share));
    }
    if let (Some(max), Some(s)) = (dp.mapq0_max, share)
        && s > max
    {
        return Ok(b.stop(Reason::HighMapq0Share));
    }
    if dp.offset > 0.0 && est <= 0.0 {
        return Ok(b.stop(Reason::BelowFloor));
    }
    if dp.ci_cv {
        let tv: Vec<f64> = tg.iter().map(|s| m[&s.name]).collect();
        let cv_c: Vec<f64> = ct.iter().map(|s| m[&s.name]).collect();
        let rel = (cv(&tv).powi(2) / tg.len() as f64 + cv(&cv_c).powi(2) / ct.len() as f64).sqrt();
        b.ci95 = Some((
            (total * (-Z95 * rel).exp() - dp.offset).max(0.0),
            total * (Z95 * rel).exp() - dp.offset,
        ));
    }
    match &row.call {
        CallRule::Psv { max_off } | CallRule::Round { max_off } => {
            let n = round_i(est);
            if (est - n as f64).abs() > *max_off {
                return Ok(b.stop(Reason::BetweenCounts));
            }
            if matches!(row.call, CallRule::Psv { .. }) {
                b.module_total = Some(n);
            } else {
                b.copies = Some(n);
                b.reading = row.reading_count(n);
            }
        }
        CallRule::Cut(cuts) => {
            let n = cut_copies(est, cuts);
            b.copies = Some(n);
            b.reading = row.reading_count(n);
        }
        CallRule::Class { low, high } => {
            let (lo, hi) = b.ci95.expect("a class row has interval cv");
            let class = if hi < *low {
                Class::Low
            } else if lo > *high {
                Class::High
            } else {
                Class::Typical
            };
            b.class = Some(class);
            b.reading = row.reading_word(class.as_str());
        }
    }
    b.status = Status::Called;
    Ok(b)
}

/// The number of cut points a value reaches (`cut:0.5,1.5` gives 0, 1 or 2).
pub fn cut_copies(value: f64, cuts: &[f64]) -> i64 {
    cuts.iter().filter(|c| value >= **c).count() as i64
}

/// Gene-base, pseudogene-base and other reads at one PSV site.
pub type SiteCounts = (u64, u64, u64);

/// A `call psv` row's block (from [`assess`]) with its sites read. With the
/// module total called, each site of [`PSV_MIN_READS`] reads or more gets
/// the copies carrying the gene's base there (total x gene reads / gene and
/// pseudogene reads) and its ratio to the median of the sites. A site is
/// low when its estimate is more than [`PSV_MARGIN`] under the larger of
/// the median's count and the usual. The CYP21A2 copies are the least count
/// over the low sites, or the median's count when none is low. Finding:
/// `deletion` when a site is low and the module total under the usual;
/// `conversion` at the low sites when the module total is kept;
/// `pseudogene_deletion` when the module total is under the usual with no
/// low site; else `none`.
pub fn psv_assess(row: &Locus, mut b: Block, counts: &HashMap<String, SiteCounts>) -> Block {
    b.sites = row
        .psv
        .iter()
        .map(|p| {
            let (g, ps, o) = counts.get(&p.name).copied().unwrap_or((0, 0, 0));
            SiteRead {
                name: p.name.clone(),
                gene_reads: g,
                pseudo_reads: ps,
                other_reads: o,
                copies: None,
                estimate: None,
                ratio: None,
            }
        })
        .collect();
    b.finding = None;
    b.conversion_sites = Vec::new();
    let (Status::Called, Some(module_total)) = (b.status, b.module_total) else {
        b.status = Status::NotAssessed;
        return b;
    };
    let total = b.total.expect("a called block has a total");
    let expected = b.expected.expect("a called block has an expected count");
    let usual = expected.div_euclid(row.targets.len() as i64);
    let mut read = Vec::new();
    for (i, s) in b.sites.iter_mut().enumerate() {
        let n = s.gene_reads + s.pseudo_reads;
        if n >= PSV_MIN_READS {
            let e = total * s.gene_reads as f64 / n as f64;
            s.estimate = Some(e);
            s.copies = Some(round_i(e));
            read.push(i);
        }
    }
    if read.is_empty() {
        b.status = Status::NotAssessed;
        b.reason = Some(Reason::PsvLowDepth);
        b.copies = None;
        return b;
    }
    let mut est: Vec<f64> = read.iter().filter_map(|&i| b.sites[i].estimate).collect();
    est.sort_by(f64::total_cmp);
    let mid = est.len() / 2;
    let median = if est.len() % 2 == 1 {
        est[mid]
    } else {
        (est[mid - 1] + est[mid]) / 2.0
    };
    for &i in &read {
        let e = b.sites[i].estimate.expect("read");
        b.sites[i].ratio = (median > 0.0).then(|| e / median);
    }
    let floor = round_i(median).max(usual) as f64;
    let conv: Vec<&SiteRead> = b
        .sites
        .iter()
        .filter(|s| s.estimate.is_some_and(|e| e < floor - PSV_MARGIN))
        .collect();
    // Only a site past the margin lowers the count, so one noisy site never
    // gives fewer copies than the finding says.
    let copies = conv
        .iter()
        .filter_map(|s| s.copies)
        .min()
        .unwrap_or_else(|| round_i(median));
    let conv_names: Vec<String> = conv.iter().map(|s| s.name.clone()).collect();
    let finding = if !conv_names.is_empty() && module_total < expected {
        Finding::Deletion
    } else if !conv_names.is_empty() {
        b.conversion_sites = conv_names.clone();
        Finding::Conversion
    } else if module_total < expected {
        Finding::PseudogeneDeletion
    } else {
        Finding::None
    };
    b.copies = Some(copies);
    b.finding = Some(finding);
    let sites = conv_names
        .iter()
        .map(|n| {
            let text = row
                .psv
                .iter()
                .find(|p| &p.name == n)
                .map_or("", |p| p.text.as_str());
            format!("{n} ({text})")
        })
        .collect::<Vec<_>>()
        .join(", ");
    b.reading = row
        .reading_word(finding.as_str())
        .replace("{sites}", &sites);
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::Catalog;

    fn catalog() -> Catalog {
        Catalog::load(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../../bin/cn_loci.tsv"
        ))
        .unwrap()
    }

    /// Depths with every segment of `row` at `x`, `over` overriding some.
    fn means(row: &Locus, x: f64, over: &[(&str, f64)]) -> Means {
        let mut m: Means = row
            .targets
            .iter()
            .chain(&row.controls)
            .map(|s| (s.name.clone(), x))
            .collect();
        for (k, v) in over {
            m.insert((*k).to_owned(), *v);
        }
        m
    }

    fn scaled(m: &Means, f: f64) -> Means {
        m.iter().map(|(k, v)| (k.clone(), v * f)).collect()
    }

    #[test]
    fn hba_counts_and_guards() {
        let cat = catalog();
        let hba = &cat.loci["HBA"];
        for (h2, h1, want) in [
            (30.0, 30.0, Some(4)),
            (15.0, 30.0, Some(3)),
            (15.0, 15.0, Some(2)),
        ] {
            let m = means(hba, 30.0, &[("HBA2", h2), ("HBA1", h1)]);
            let b = assess(hba, &m, Some(&scaled(&m, 0.75)), Sex::Unknown).unwrap();
            assert_eq!((b.status, b.copies), (Status::Called, want));
        }
        let m = means(hba, 30.0, &[("HBA2", 22.5)]);
        let b = assess(hba, &m, Some(&scaled(&m, 0.75)), Sex::Unknown).unwrap();
        assert_eq!(b.reason, Some(Reason::BetweenCounts));
        let m = means(hba, 30.0, &[]);
        let b = assess(hba, &m, Some(&scaled(&m, 0.99)), Sex::Unknown).unwrap();
        assert_eq!(b.reason, Some(Reason::LowMapq0Share));
        let b = assess(hba, &Means::new(), Some(&Means::new()), Sex::Unknown).unwrap();
        assert_eq!(b.reason, Some(Reason::NoData));
    }

    #[test]
    fn a_control_with_no_read_is_a_broken_control() {
        let cat = catalog();
        let hba = &cat.loci["HBA"];
        let m = means(hba, 30.0, &[("flank_left", 0.0)]);
        let e = assess(hba, &m, Some(&m), Sex::Unknown).unwrap_err();
        assert!(
            e.to_string()
                .starts_with("BROKEN_CONTROL: control segment flank_left of HBA"),
            "{e}"
        );
    }

    #[test]
    fn ties_round_to_even() {
        // An HBA estimate of exactly 2.5 lies 0.5 from both counts:
        // between_counts whatever the rounding; with round:0.5 it would call
        // 2 (half to even), never 3.
        assert_eq!(round_i(2.5), 2);
        assert_eq!(round_i(3.5), 4);
        assert_eq!(round_i(0.5), 0);
    }

    #[test]
    fn chrx_and_chry_follow_the_sex() {
        let cat = catalog();
        let az = cat
            .loci
            .values()
            .find(|l| l.chrom == "chrY" && l.is_interval_row())
            .unwrap();
        let m = means(az, 15.0, &[]);
        let b = assess(az, &m, None, Sex::Female).unwrap();
        assert_eq!(
            (b.status, b.reason),
            (Status::NotApplicable, Some(Reason::NoY))
        );
        let b = assess(az, &m, None, Sex::Unknown).unwrap();
        assert_eq!(b.reason, Some(Reason::SexUnknown));
        let b = assess(az, &m, None, Sex::Male).unwrap();
        let units = if az.target_depth == "sum" {
            az.targets.len() as i64
        } else {
            1
        };
        assert_eq!(b.expected, Some(units));
    }

    #[test]
    fn lpa_class_and_interval() {
        let cat = catalog();
        let lpa = &cat.loci["LPA_KIV2"];
        // N + 2 = 14 x M / D: M = 30 / 14 x (N + 2) at D = 30.
        let m = means(lpa, 30.0, &[]);
        let kiv: Vec<(String, f64)> = lpa
            .targets
            .iter()
            .map(|s| (s.name.clone(), 30.0 / 14.0 * 22.0))
            .collect();
        let mut m2 = m.clone();
        m2.extend(kiv);
        let b = assess(lpa, &m2, Some(&scaled(&m2, 0.1)), Sex::Unknown).unwrap();
        assert_eq!(b.class, Some(Class::Low));
        let (lo, hi) = b.ci95.unwrap();
        assert!((b.estimate.unwrap() - 20.0).abs() < 1e-9 && lo <= 20.0 && hi >= 20.0);
    }
}
