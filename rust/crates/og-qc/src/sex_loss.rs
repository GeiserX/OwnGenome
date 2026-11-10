//! `og qc sex-loss`: the share of blood cells that lost the Y (a man) or
//! one X (a woman), added as `mosaic_*` rows to step 33's key/value table.
//!
//! ```text
//! og qc sex-loss --table S_sample_qc.tsv [--regions S.regions.bed.gz] [--vcf S.vcf.gz]
//! ```
//!
//! The sex is the one somalier read from the reads (`inferred_sex` in the
//! table), else the declared one; with neither, nothing is estimated.
//!
//! **Depth.** mosdepth's 500 bp windows are grouped into 10 kb blocks; a
//! block is the median of its windows with reads, and needs 10 of them. The
//! ratio is the median block of the sex chromosome over the median block of
//! chr1-22. Only the single-copy part of the sex chromosome counts
//! ([`SINGLE_COPY_Y`], [`SINGLE_COPY_X`]): a whole-chromosome mean would
//! count chrY's gaps, its repeats and the copies of chrX it holds. One Y in
//! every cell reads 0.5, two X read 1.0. With a share `f` of the cells
//! missing a copy the ratio drops by `f / k` of its expected value (`k` the
//! normal number of copies), so `f = k (1 - ratio / expected)`.
//!
//! **Allele fractions (women only).** At heterozygous chrX SNVs, one X lost
//! in a share `f` of the cells moves the allele fraction from 0.5 to
//! `0.5 ± d`, `d = f / (2 (2 - f))`, so `f = 4d / (1 + 2d)`. Per site,
//! `u = (alt - ref)² / n` has mean `φ + 4 (n - 1) d²`, `φ` the library's
//! spread at that depth, which chr1-22 measure:
//! `d² = (mean u on chrX - mean u on chr1-22) / (4 (mean n - 1))`.
//!
//! **Uncertainty**, as a 95% range: each median's standard error
//! (`1.2533 × 1.4826 × MAD / √blocks`) and each mean's (from the sample
//! variance). Depth also drifts along the genome with GC content and
//! mappability, and neighbouring blocks share that drift, so the depth range
//! adds the spread (`1.4826 × MAD`) of the medians of 200 kb stretches of
//! chr1-22, over the square root of the number of stretches the sex
//! chromosome's blocks fill. A loss counts as detected only when the low end
//! of its range is above [`FLOOR`]: at 30x the sex-chromosome depth varies
//! by a few percent between libraries.
//!
//! The arithmetic is CPython's, step by step: the medians, `fmean` and the
//! exact `variance` of `pycompat::stats`, and every quotient correctly
//! rounded, so the rows match the Python this replaced to the last digit.

use og_io::{Lines, Utf8};
use pycompat::{fmt, json::Value, num, stats, text, value};
use std::collections::HashMap;
use std::io::{self, Write};

/// Below this share of the cells nothing can be told at 30x.
pub const FLOOR: f64 = 0.10;
/// bp per block of mosdepth windows.
const BLOCK: i128 = 10_000;
/// bp per stretch, for the drift of depth along the genome.
const STRETCH: i128 = 200_000;
/// Blocks per stretch.
const PER: usize = 20;
/// Windows with reads a block needs (a block holds 20 at `--by 500`).
const MIN_WINDOWS: usize = 10;
/// Blocks the sex chromosome and the autosomes each need.
const MIN_BLOCKS: usize = 10;
/// Stretches of chr1-22 (with `MIN_BLOCKS` blocks each) the drift needs.
const MIN_STRETCHES: usize = 5;
/// Reads a heterozygous site needs.
const MIN_DP: i128 = 10;
/// Heterozygous chrX sites the allele-fraction estimate needs.
const MIN_X_HETS: usize = 100;
/// Heterozygous autosomal sites it needs for the library's spread.
const MIN_AUTO_HETS: usize = 1000;
const Z95: f64 = 1.96;

/// chrY's single-copy part on GRCh38, 0-based half-open: GIAB genome
/// stratifications v3.3 (GRCh38@all/XY), `chrY_nonPAR` minus `_XTR` and
/// `_ampliconic`, minus the centromere (10,316,944-10,544,039), pieces of
/// 200 kb or more.
pub const SINGLE_COPY_Y: &[(i128, i128)] = &[
    (2781479, 3050044),
    (6748713, 7574481),
    (11037032, 11592902),
    (11674123, 13984473),
    (14058733, 15875093),
    (15905214, 16159795),
    (16425965, 17456266),
    (18870334, 20207793),
    (20257793, 21305953),
    (26415435, 26673214),
];

/// chrX's single-copy part, the same way (centromere 58,100,000-61,000,000).
pub const SINGLE_COPY_X: &[(i128, i128)] = &[
    (2781479, 37099262),
    (37285837, 48343310),
    (48434604, 49157254),
    (49528394, 50228964),
    (50278964, 52032464),
    (52223402, 52489200),
    (52520537, 55437684),
    (55547739, 58100000),
    (61000000, 62412542),
    (62462542, 63185573),
    (63193371, 71674267),
    (71835832, 72721314),
    (73105236, 89140830),
    (93428068, 101563740),
    (101648990, 102180805),
    (102519463, 103940531),
    (104117650, 114281198),
    (114331198, 115738949),
    (115838949, 116557779),
    (116595566, 120879381),
    (120929381, 139870917),
    (139871081, 144425606),
    (144475606, 155701383),
];

/// The rows this tool owns, in the order it writes them.
pub const KEYS: [&str; 16] = [
    "mosaic_sex",
    "mosaic_sex_from",
    "mosaic_chrom",
    "mosaic_depth_ratio",
    "mosaic_depth_ratio_expected",
    "mosaic_depth_blocks",
    "mosaic_loss_depth",
    "mosaic_loss_depth_low",
    "mosaic_loss_depth_high",
    "mosaic_x_het_sites",
    "mosaic_loss_af",
    "mosaic_loss_af_low",
    "mosaic_loss_af_high",
    "mosaic_floor",
    "mosaic_call",
    "mosaic_reason",
];

/// The sex whose sex chromosome is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
}

impl Sex {
    fn parse(s: Option<&str>) -> Option<Sex> {
        match s {
            Some("male") => Some(Sex::Male),
            Some("female") => Some(Sex::Female),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Sex::Male => "male",
            Sex::Female => "female",
        }
    }

    /// The chromosome read, without `chr`.
    fn chrom(self) -> &'static str {
        match self {
            Sex::Male => "Y",
            Sex::Female => "X",
        }
    }

    /// Its depth ratio to chr1-22 with no loss.
    fn expected(self) -> f64 {
        match self {
            Sex::Male => 0.5,
            Sex::Female => 1.0,
        }
    }

    /// Its normal number of copies.
    fn copies(self) -> f64 {
        match self {
            Sex::Male => 1.0,
            Sex::Female => 2.0,
        }
    }

    fn pieces(self) -> &'static [(i128, i128)] {
        match self {
            Sex::Male => SINGLE_COPY_Y,
            Sex::Female => SINGLE_COPY_X,
        }
    }
}

/// True when `[start, end)` lies inside one piece of `pieces`.
fn single_copy(pieces: &[(i128, i128)], start: i128, end: i128) -> bool {
    let i = pieces.partition_point(|&(s, _)| s <= start);
    i > 0 && end <= pieces[i - 1].1
}

/// The chromosome's name without a leading `chr`.
fn bare(chrom: &str) -> &str {
    chrom.strip_prefix("chr").unwrap_or(chrom)
}

/// `1` to `22`, written as Python's `str(int)` writes them.
fn is_autosome(name: &str) -> bool {
    let b = name.as_bytes();
    match b {
        [d] => (b'1'..=b'9').contains(d),
        [a, d] => (*a == b'1' && d.is_ascii_digit()) || (*a == b'2' && (b'0'..=b'2').contains(d)),
        _ => false,
    }
}

/// Why an input could not be read, as the reason row names it: a Python
/// exception class and its message, the same text the Python wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadError {
    class: &'static str,
    message: String,
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.class, self.message)
    }
}

/// `repr(s)` of a str.
fn repr(s: &str) -> String {
    value::py_repr(&Value::Str(s.to_owned())).unwrap_or_else(|_| format!("'{s}'"))
}

impl ReadError {
    fn value(message: String) -> Self {
        ReadError {
            class: "ValueError",
            message,
        }
    }

    /// An error opening or reading `path`, `line` the line being read.
    fn io(e: &io::Error, path: &str, line: u64) -> Self {
        if let Some(code) = e.raw_os_error() {
            let (class, text) = match code {
                2 => ("FileNotFoundError", "No such file or directory".to_owned()),
                13 => ("PermissionError", "Permission denied".to_owned()),
                21 => ("IsADirectoryError", "Is a directory".to_owned()),
                _ => {
                    let s = e.to_string();
                    let text = s.split(" (os error").next().unwrap_or(&s).to_owned();
                    ("OSError", text)
                }
            };
            return ReadError {
                class,
                message: format!("[Errno {code}] {text}: {}", repr(path)),
            };
        }
        let utf8 = e
            .get_ref()
            .is_some_and(|inner| inner.is::<std::str::Utf8Error>());
        match e.kind() {
            io::ErrorKind::UnexpectedEof => ReadError {
                class: "EOFError",
                message: "Compressed file ended before the end-of-stream marker was reached"
                    .to_owned(),
            },
            _ if utf8 => ReadError {
                class: "UnicodeDecodeError",
                message: format!("line {line} is not UTF-8"),
            },
            // Python's zlib.error; the text is the decoder's own.
            _ => ReadError {
                class: "error",
                message: format!("line {line}: {e}"),
            },
        }
    }
}

/// `int(s)`. Plain digits, as mosdepth writes them, skip the general parser.
fn int(s: &str) -> Result<i128, ReadError> {
    if (1..=18).contains(&s.len()) && s.bytes().all(|c| c.is_ascii_digit()) {
        return Ok(s.bytes().fold(0, |v, c| v * 10 + i128::from(c - b'0')));
    }
    num::parse_int(s).map_err(|e| match e {
        num::ParseError::Invalid => ReadError::value(format!(
            "invalid literal for int() with base 10: {}",
            repr(s)
        )),
        num::ParseError::TooLarge => {
            ReadError::value(format!("{} does not fit 128 bits", repr(s.trim())))
        }
    })
}

/// `float(s)`. `12.34` with or without a newline, as mosdepth writes it,
/// skips the general parser: Rust's own parse is correctly rounded too.
fn float(s: &str) -> Result<f64, ReadError> {
    let t = s.strip_suffix('\n').unwrap_or(s);
    let plain = |p: &str| !p.is_empty() && p.bytes().all(|c| c.is_ascii_digit());
    let fast = match t.split_once('.') {
        Some((a, b)) => plain(a) && plain(b),
        None => plain(t),
    };
    if fast && let Ok(v) = t.parse::<f64>() {
        return Ok(v);
    }
    num::parse_float(s)
        .map_err(|_| ReadError::value(format!("could not convert string to float: {}", repr(s))))
}

/// Opens `path` (plain, gzip or BGZF, as mosdepth and bcftools write it,
/// inflated on every core) as Python's text mode reads it.
fn lines(path: &str) -> Result<Lines<Box<dyn io::BufRead + Send>>, ReadError> {
    let r =
        og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(|e| ReadError::io(&e, path, 0))?;
    Ok(Lines::new(r, Utf8::Strict))
}

/// `statistics.median`: the data here are never empty and always finite.
fn median(xs: impl IntoIterator<Item = f64>) -> f64 {
    stats::median(xs).unwrap_or(f64::NAN)
}

/// The median and its standard error (from the MAD).
fn median_se(xs: &[f64]) -> (f64, f64) {
    let m = median(xs.iter().copied());
    let mad = median(xs.iter().map(|x| (x - m).abs()));
    (m, 1.2533 * 1.4826 * mad / (xs.len() as f64).sqrt())
}

/// The block medians of a mosdepth regions file.
#[derive(Debug, Default)]
pub struct Blocks {
    /// chr1-22, in file order.
    pub auto: Vec<f64>,
    /// chr1-22 grouped by their 200 kb stretch.
    pub stretches: Vec<Vec<f64>>,
    /// The single-copy part of the sex chromosome.
    pub sex: Vec<f64>,
}

/// The current block while the windows stream by.
struct Open {
    is_sex: bool,
    contig: String,
    block: i128,
    depths: Vec<f64>,
}

impl Blocks {
    /// The median of `open`'s windows, when it has enough of them.
    fn close(&mut self, open: &Open, index: &mut HashMap<(String, i128), usize>) {
        if open.depths.len() < MIN_WINDOWS {
            return;
        }
        let m = median(open.depths.iter().copied());
        if open.is_sex {
            self.sex.push(m);
            return;
        }
        self.auto.push(m);
        let key = (
            open.contig.clone(),
            (open.block * BLOCK).div_euclid(STRETCH),
        );
        let next = self.stretches.len();
        let i = *index.entry(key).or_insert(next);
        if i == next {
            self.stretches.push(Vec::new());
        }
        self.stretches[i].push(m);
    }
}

/// Reads mosdepth's windows: block medians on chr1-22 and on the
/// single-copy part of `sex`'s chromosome. Windows with no reads are gaps,
/// or outside what a sliced or targeted BAM holds, and do not count.
pub fn depth_blocks(path: &str, sex: Sex) -> Result<Blocks, ReadError> {
    let mut lines = lines(path)?;
    let mut blocks = Blocks::default();
    let mut index = HashMap::new();
    let mut open: Option<Open> = None;
    let mut line = String::new();
    let mut n = 0u64;
    loop {
        match lines.read_line(&mut line) {
            Ok(true) => n += 1,
            Ok(false) => break,
            Err(e) => return Err(ReadError::io(&e, path, n + 1)),
        }
        let mut cols = line.splitn(5, '\t');
        let (Some(c), Some(s), Some(e), Some(d)) =
            (cols.next(), cols.next(), cols.next(), cols.next())
        else {
            let got = line.splitn(5, '\t').count();
            return Err(ReadError::value(format!(
                "not enough values to unpack (expected 4, got {got})"
            )));
        };
        let name = bare(c);
        let is_sex = if is_autosome(name) {
            false
        } else if name == sex.chrom() {
            let (start, end) = (int(s)?, int(e)?);
            if !single_copy(sex.pieces(), start, end) {
                continue;
            }
            true
        } else {
            continue;
        };
        let dp = float(d)?;
        if dp <= 0.0 {
            continue;
        }
        if !dp.is_finite() {
            return Err(ReadError::value(format!(
                "line {n}: a depth of {}, not a finite number",
                fmt::repr(dp)
            )));
        }
        let block = int(s)?.div_euclid(BLOCK);
        let same = open
            .as_ref()
            .is_some_and(|o| o.is_sex == is_sex && o.block == block && o.contig == c);
        if !same {
            if let Some(o) = open.as_mut() {
                blocks.close(o, &mut index);
                o.depths.clear();
                o.contig.clear();
                o.contig.push_str(c);
                o.is_sex = is_sex;
                o.block = block;
            } else {
                open = Some(Open {
                    is_sex,
                    contig: c.to_owned(),
                    block,
                    depths: Vec::with_capacity(32),
                });
            }
        }
        if let Some(o) = open.as_mut() {
            o.depths.push(dp);
        }
    }
    if let Some(o) = open.as_ref() {
        blocks.close(o, &mut index);
    }
    Ok(blocks)
}

/// The depth estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Depth {
    pub ratio: f64,
    pub blocks: usize,
    pub loss: f64,
    pub low: f64,
    pub high: f64,
}

/// The depth estimate, or why there is none.
pub fn depth_estimate(path: &str, sex: Sex) -> Result<Result<Depth, String>, ReadError> {
    let b = depth_blocks(path, sex)?;
    let (auto, sx) = (&b.auto, &b.sex);
    if auto.len() < MIN_BLOCKS || sx.len() < MIN_BLOCKS {
        return Ok(Err(format!(
            "too few 10 kb blocks with reads: {} on chr{}'s single-copy part and {} on chr1-22 \
             (each needs {MIN_BLOCKS}; a targeted or sliced BAM may have none)",
            sx.len(),
            sex.chrom(),
            auto.len()
        )));
    }
    let (ma, sea) = median_se(auto);
    let (ms, ses) = median_se(sx);
    let levels: Vec<f64> = b
        .stretches
        .iter()
        .filter(|v| v.len() >= MIN_BLOCKS)
        .map(|v| median(v.iter().copied()) / ma)
        .collect();
    if levels.len() < MIN_STRETCHES {
        return Ok(Err(format!(
            "{} stretches of 200 kb with {MIN_BLOCKS} or more blocks on chr1-22 \
             (needs {MIN_STRETCHES} to measure how depth drifts along the genome)",
            levels.len()
        )));
    }
    let mid = median(levels.iter().copied());
    let drift = 1.4826 * median(levels.iter().map(|x| (x - mid).abs()));
    let per = PER as f64;
    let r = ms / ma;
    let sq = |x: f64| x * x;
    let se_r = if ms > 0.0 {
        r * (sq(ses / ms)
            + sq(sea / ma)
            + sq(drift) * (per / sx.len().max(PER) as f64 + per / auto.len().max(PER) as f64))
            .sqrt()
    } else {
        0.0
    };
    let f = sex.copies() * (1.0 - r / sex.expected());
    let se_f = sex.copies() * se_r / sex.expected();
    Ok(Ok(Depth {
        ratio: r,
        blocks: sx.len(),
        loss: f,
        low: f - Z95 * se_f,
        high: f + Z95 * se_f,
    }))
}

/// `num / den` correctly rounded, as Python's int true division gives it.
/// `den` is at most 2^32 here.
fn ratio(num: u128, den: u128) -> f64 {
    if num < 1 << 53 {
        // Both exact as floats, so one division rounds once.
        return num as f64 / den as f64;
    }
    // Scale so the integer quotient has 55 bits or more, round it to 53
    // with the remainder as a sticky bit, then scale back (exact).
    let shift = (den.ilog2() + 56).saturating_sub(num.ilog2());
    let n = num << shift;
    let (q, rem) = (n / den, n % den);
    let drop = q.ilog2() + 1 - 53;
    let mut m = q >> drop;
    let rest = q & ((1 << drop) - 1);
    let half = 1 << (drop - 1);
    if rest > half || (rest == half && (rem != 0 || m & 1 == 1)) {
        m += 1;
    }
    m as f64 * 2f64.powi(drop as i32) / 2f64.powi(shift as i32)
}

/// The PASS heterozygous SNVs of chr1-22 and of chrX's single-copy part:
/// each site's read count `n` and `u = (alt - ref)² / n`.
#[derive(Debug, Default)]
pub struct Hets {
    pub auto_n: Vec<i64>,
    pub auto_u: Vec<f64>,
    pub x_n: Vec<i64>,
    pub x_u: Vec<f64>,
}

/// One site's `(n, u)` from its first sample's `GT` and `AD`, or `None`
/// when it is not a heterozygous site with enough reads.
fn het(format: &str, sample: &str) -> Option<(i64, f64)> {
    let keys: Vec<&str> = format.split(':').collect();
    let vals: Vec<&str> = sample.split(':').collect();
    let gt = keys.iter().position(|&k| k == "GT")?;
    let ad = keys.iter().position(|&k| k == "AD")?;
    if vals.len() < keys.len() || !matches!(vals[gt], "0/1" | "0|1" | "1/0" | "1|0") {
        return None;
    }
    let mut ad = vals[ad].split(',');
    let r = num::parse_int(ad.next()?).ok()?;
    let a = num::parse_int(ad.next()?).ok()?;
    let n = r.checked_add(a)?;
    // A read count past 32 bits is no read count: the Python stopped with
    // an OverflowError there, og skips the site.
    if n < MIN_DP || n > i128::from(u32::MAX) {
        return None;
    }
    let diff = a.checked_sub(r)?.unsigned_abs();
    if diff >> 64 != 0 {
        return None;
    }
    Some((n as i64, ratio(diff * diff, n as u128)))
}

/// Reads the heterozygous SNVs of a VCF.
pub fn het_sites(path: &str) -> Result<Hets, ReadError> {
    let mut lines = lines(path)?;
    let mut h = Hets::default();
    let mut line = String::new();
    let mut n = 0u64;
    loop {
        match lines.read_line(&mut line) {
            Ok(true) => n += 1,
            Ok(false) => break,
            Err(e) => return Err(ReadError::io(&e, path, n + 1)),
        }
        if line.starts_with('#') {
            continue;
        }
        let mut p = [""; 11];
        let mut k = 0;
        for (slot, col) in p.iter_mut().zip(line.splitn(11, '\t')) {
            *slot = col;
            k += 1;
        }
        if k < 10 {
            continue;
        }
        let name = bare(p[0]);
        let to_x = if is_autosome(name) {
            false
        } else if name == "X" {
            let pos = int(p[1])?;
            if !single_copy(SINGLE_COPY_X, pos - 1, pos) {
                continue;
            }
            true
        } else {
            continue;
        };
        if !matches!(p[6], "PASS" | ".")
            || p[3].chars().count() != 1
            || p[4].chars().count() != 1
            || matches!(p[4], "." | "*")
        {
            continue;
        }
        let Some((dp, u)) = het(p[8], text::rstrip(p[9], Some("\n"))) else {
            continue;
        };
        let (ns, us) = if to_x {
            (&mut h.x_n, &mut h.x_u)
        } else {
            (&mut h.auto_n, &mut h.auto_u)
        };
        ns.push(dp);
        us.push(u);
    }
    Ok(h)
}

/// The allele-fraction estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Af {
    pub sites: usize,
    pub loss: f64,
    pub low: f64,
    pub high: f64,
}

/// `max(v, 0.0)` as Python takes it: `v` unless `0.0 > v` (so NaN and -0.0 stay).
fn py_max0(v: f64) -> f64 {
    if 0.0 > v { 0.0 } else { v }
}

/// The share of cells `f` that the allele-fraction shift `d²` reads as.
fn af_loss(d2: f64) -> f64 {
    let d = py_max0(d2).sqrt();
    4.0 * d / (1.0 + 2.0 * d)
}

/// The allele-fraction estimate, or why there is none.
pub fn af_estimate(path: &str) -> Result<Result<Af, String>, ReadError> {
    let h = het_sites(path)?;
    if h.auto_n.len() < MIN_AUTO_HETS {
        return Ok(Err(format!(
            "{} heterozygous SNVs on chr1-22 with {MIN_DP} or more reads (needs {MIN_AUTO_HETS})",
            h.auto_n.len()
        )));
    }
    // Sites far deeper than the genome are collapsed copies; their allele
    // fractions are not one locus'.
    let cap = stats::median_int(h.auto_n.iter().copied()).unwrap_or(stats::Median::Float(0.0));
    let keep = |n: i64| match cap {
        stats::Median::Int(m) => i128::from(n) <= 2 * i128::from(m),
        stats::Median::Float(m) => (n as f64) <= 2.0 * m,
    };
    let a: Vec<f64> = h
        .auto_n
        .iter()
        .zip(&h.auto_u)
        .filter(|&(&n, _)| keep(n))
        .map(|(_, &u)| u)
        .collect();
    let (xn, xs): (Vec<i64>, Vec<f64>) = h
        .x_n
        .iter()
        .zip(&h.x_u)
        .filter(|&(&n, _)| keep(n))
        .map(|(&n, &u)| (n, u))
        .unzip();
    if xs.len() < MIN_X_HETS {
        return Ok(Err(format!(
            "{} heterozygous SNVs on chrX's single-copy part (needs {MIN_X_HETS})",
            xs.len()
        )));
    }
    // Never empty, never fewer than two, never overflowing here.
    let fmean = |v: &[f64]| stats::fmean(v.iter().copied()).unwrap_or(f64::NAN);
    // RED CONTROL, reverted in the next commit: a naive two-pass variance.
    let var = |v: &[f64]| {
        let mu = v.iter().sum::<f64>() / v.len() as f64;
        v.iter().map(|x| (x - mu) * (x - mu)).sum::<f64>() / (v.len() - 1) as f64
    };
    let mean_n = stats::fmean(xn.iter().map(|&n| n as f64)).unwrap_or(f64::NAN);
    let scale = 4.0 * (mean_n - 1.0);
    let d2 = (fmean(&xs) - fmean(&a)) / scale;
    let se = (var(&xs) / xs.len() as f64 + var(&a) / a.len() as f64).sqrt() / scale;
    Ok(Ok(Af {
        sites: xs.len(),
        loss: af_loss(d2),
        low: af_loss(d2 - Z95 * se),
        high: af_loss(d2 + Z95 * se),
    }))
}

fn fmt3(v: f64) -> String {
    fmt::fixed(v, 3)
}

/// `os.path.basename`.
fn basename(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

/// The mosaic rows for one sample, in [`KEYS`] order, from its table's
/// other rows and the two inputs (`None` or empty when there is none).
pub fn estimate(
    table: &[(String, String)],
    regions: Option<&str>,
    vcf: Option<&str>,
) -> Vec<(&'static str, String)> {
    let mut row: Vec<(&'static str, String)> = KEYS.iter().map(|&k| (k, String::new())).collect();
    let mut set = |key: &str, v: String| {
        if let Some(slot) = row.iter_mut().find(|(k, _)| *k == key) {
            slot.1 = v;
        }
    };
    // The last row of a key wins, as dict() keeps it.
    let get = |key: &str| {
        table
            .iter()
            .rev()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    };
    let (sex, source) = match Sex::parse(get("inferred_sex")) {
        Some(s) => (Some(s), "reads"),
        None => (Sex::parse(get("declared_sex")), "declared"),
    };
    set("mosaic_floor", fmt3(FLOOR));
    let Some(sex) = sex else {
        set("mosaic_call", "not_estimated".to_owned());
        set(
            "mosaic_reason",
            "no sex: none declared and somalier could not tell it from the reads".to_owned(),
        );
        return row;
    };
    let chrom = sex.chrom();
    set("mosaic_sex", sex.name().to_owned());
    set("mosaic_sex_from", source.to_owned());
    set("mosaic_chrom", format!("chr{chrom}"));
    set("mosaic_depth_ratio_expected", fmt::fixed(sex.expected(), 1));
    let mut why: Vec<String> = Vec::new();
    let (mut dep_lo, mut dep_hi, mut af_lo) = (None, None, None);
    match regions.filter(|p| !p.is_empty()) {
        None => why.push(
            "no mosdepth windows (step 16b's <sample>.regions.bed.gz): depth not read".to_owned(),
        ),
        Some(path) => match depth_estimate(path, sex) {
            Ok(Ok(d)) => {
                set("mosaic_depth_ratio", fmt3(d.ratio));
                set("mosaic_depth_blocks", d.blocks.to_string());
                set("mosaic_loss_depth", fmt3(d.loss));
                set("mosaic_loss_depth_low", fmt3(d.low));
                set("mosaic_loss_depth_high", fmt3(d.high));
                (dep_lo, dep_hi) = (Some(d.low), Some(d.high));
            }
            Ok(Err(reason)) => why.push(format!("depth not estimated: {reason}")),
            Err(e) => why.push(format!(
                "depth not estimated: could not read {}: {e}",
                basename(path)
            )),
        },
    }
    if sex == Sex::Female {
        match vcf.filter(|p| !p.is_empty()) {
            None => why.push("no VCF: chrX allele fractions not read".to_owned()),
            Some(path) => match af_estimate(path) {
                Ok(Ok(a)) => {
                    set("mosaic_x_het_sites", a.sites.to_string());
                    set("mosaic_loss_af", fmt3(a.loss));
                    set("mosaic_loss_af_low", fmt3(a.low));
                    set("mosaic_loss_af_high", fmt3(a.high));
                    af_lo = Some(a.low);
                }
                Ok(Err(reason)) => why.push(format!("allele fractions not estimated: {reason}")),
                Err(e) => why.push(format!(
                    "allele fractions not estimated: could not read {}: {e}",
                    basename(path)
                )),
            },
        }
    }
    let what = match sex {
        Sex::Male => "the Y",
        Sex::Female => "one X",
    };
    let floor = fmt::percent(FLOOR, 0);
    let call = if dep_lo.is_none() && af_lo.is_none() {
        "not_estimated"
    } else if dep_hi.is_some_and(|h| h < -FLOOR) {
        // An extra copy also shifts chrX allele fractions (1/3 and 2/3 in 47,XXX).
        let copies = match sex {
            Sex::Male => "one copy",
            Sex::Female => "two copies",
        };
        why.insert(
            0,
            format!(
                "more chr{chrom} depth than {copies}: an extra copy (47,XYY or 47,XXX) or a gain, not a loss"
            ),
        );
        "gain"
    } else if dep_lo.is_some_and(|l| l > FLOOR) {
        why.insert(
            0,
            format!(
                "the low end of the 95% range is above {floor}: a share of the cells lost {what}"
            ),
        );
        "detected"
    } else if af_lo.is_some_and(|l| l > FLOOR) {
        why.insert(
            0,
            format!(
                "chrX allele fractions are shifted beyond {floor} but its depth is not: in a share of \
                 the cells one X was lost and the other copied (copy-neutral), or the depth cannot tell"
            ),
        );
        "detected"
    } else {
        why.insert(
            0,
            format!(
                "the 95% range does not put a loss of {what} above {floor}; \
                 below that nothing can be told at 30x"
            ),
        );
        "not_detected"
    };
    set("mosaic_call", call.to_owned());
    set("mosaic_reason", why.join("; "));
    row
}

/// Step 33's key/value table: its rows after the `key\tvalue` header, a row
/// without a tab read as a key with an empty value. Blank lines are skipped.
pub fn read_table(path: &str) -> Result<Vec<(String, String)>, String> {
    let file = std::fs::File::open(path).map_err(|e| ReadError::io(&e, path, 0).message)?;
    let mut lines = Lines::new(io::BufReader::new(file), Utf8::Strict);
    let mut rows: Vec<(String, Option<String>)> = Vec::new();
    let mut line = String::new();
    loop {
        match lines.read_line(&mut line) {
            Ok(true) => {}
            Ok(false) => break,
            Err(e) => {
                return Err(ReadError::io(&e, path, rows.len() as u64 + 1).message);
            }
        }
        if text::strip(&line, None).is_empty() {
            continue;
        }
        let l = text::rstrip(&line, Some("\n"));
        rows.push(match l.split_once('\t') {
            Some((k, v)) => (k.to_owned(), Some(v.to_owned())),
            None => (l.to_owned(), None),
        });
    }
    match rows.first() {
        Some((k, Some(v))) if k == "key" && v == "value" => {}
        _ => return Err(format!("{path} is not a key/value table")),
    }
    Ok(rows
        .into_iter()
        .skip(1)
        .map(|(k, v)| (k, v.unwrap_or_default()))
        .collect())
}

const USAGE: &str = "usage: og qc sex-loss [-h] --table TABLE [--regions REGIONS] [--vcf VCF]\n";

const HELP: &str = "
Mosaic loss of chromosome Y (men) or of one X (women) in the cells the DNA
came from, estimated from mosdepth's windows and from the VCF's allele
fractions. Adds the mosaic_* rows to step 33's table.

options:
  -h, --help         show this help message and exit
  --table TABLE      the step 33 table; the rows are added to it
  --regions REGIONS  mosdepth's <sample>.regions.bed.gz (step 16b, --by 500)
  --vcf VCF          the sample's VCF (GT and AD), for the chrX allele fractions of a woman
";

fn usage_error(msg: &str) -> u8 {
    eprint!("{USAGE}");
    eprintln!("og qc sex-loss: error: {msg}");
    2
}

/// `og qc sex-loss ARGS`: the exit code. 1 when the table cannot be read or
/// written; a missing or unreadable regions file or VCF is a reason row.
pub fn cli(args: &[String]) -> u8 {
    let (mut table, mut regions, mut vcf) = (None, None, None);
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        if arg == "-h" || arg == "--help" {
            print!("{USAGE}{HELP}");
            return 0;
        }
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_owned())),
            _ => (arg.as_str(), None),
        };
        let slot = match flag {
            "--table" => &mut table,
            "--regions" => &mut regions,
            "--vcf" => &mut vcf,
            _ => return usage_error(&format!("unrecognized arguments: {arg}")),
        };
        let value = match inline {
            Some(v) => v,
            None => match it.next() {
                Some(v) if !(v.starts_with('-') && v.len() > 1) => v.clone(),
                _ => return usage_error(&format!("argument {flag}: expected one argument")),
            },
        };
        *slot = Some(value);
    }
    let Some(table) = table else {
        return usage_error("the following arguments are required: --table");
    };
    let rows = match read_table(&table) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("ERROR: {e}");
            return 1;
        }
    };
    // A rerun replaces the rows of the last one.
    let kept: Vec<(String, String)> = rows
        .into_iter()
        .filter(|(k, _)| !k.starts_with("mosaic_"))
        .collect();
    let mosaic = estimate(&kept, regions.as_deref(), vcf.as_deref());
    let write = || -> io::Result<()> {
        let mut out = og_io::AtomicFile::create(&table)?;
        out.write_all(b"key\tvalue\n")?;
        for (k, v) in &kept {
            writeln!(out, "{k}\t{v}")?;
        }
        for (k, v) in &mosaic {
            writeln!(out, "{k}\t{v}")?;
        }
        out.commit()
    };
    if let Err(e) = write() {
        eprintln!("ERROR: could not write {table}: {e}");
        return 1;
    }
    let mut stdout = io::stdout().lock();
    for (k, v) in &mosaic {
        if writeln!(stdout, "  {k}: {v}").is_err() {
            return 1;
        }
    }
    0
}

#[cfg(test)]
mod tests;
