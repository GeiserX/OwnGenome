//! `og y private-snps` (step 37): the chrY SNVs of a sample that no public
//! SNP table names.
//!
//! ```text
//! og y private-snps --vcf VCF --tree-dir DIR --sample S --out-prefix P
//!     [--callable-mask FASTA] [--min-depth N] [--min-gq N]
//!     [--min-fraction F] [--max-depth-factor F]
//! ```
//!
//! A haplogroup call places a paternal line on the deepest branch the public
//! trees already have. The variants that would place it further down are
//! the ones on no public list yet. A record is kept when it passes every
//! gate, in this order, and each gate counts what it drops:
//!
//! 1. chrY (or `Y`) outside the pseudoautosomal regions;
//! 2. inside the 1000 Genomes strict accessibility mask of chrY (base `P`);
//! 3. FILTER `PASS`, `.` or empty;
//! 4. the genotype is the ALT allele alone (`1`, `1/1`, `1|1`);
//! 5. depth (AD summed, else DP) at least `--min-depth` (10);
//! 6. GQ at least `--min-gq` (20);
//! 7. ALT fraction (AD, else VAF) at least `--min-fraction` (0.9);
//! 8. depth at most `--max-depth-factor` (2.0) times the median depth of
//!    the calls that passed 1 to 7 (0 turns this gate off);
//! 9. not at any position of the two public SNP tables of the tree folder.
//!
//! Records with more than one REF or ALT base never reach the gates.
//! `--callable-mask` defaults to `DIR/callable_mask_chrY.fasta.gz`.
//!
//! How it is built: the mask is one bit per base, read in one pass over the
//! FASTA (57 Mb of chrY in about 7 MB); the public positions are a sorted
//! list searched by bisection; the VCF is read one line at a time, and only
//! the calls that pass gates 1 to 7 are kept, for the median.
//!
//! What changes on purpose: a mask FASTA with a second record is refused
//! (the Python joined the records into one sequence, which shifts every
//! position after the first); a gzip input is recognised by its bytes, not
//! by a `.gz` name; the output is written whole or not at all; digits are
//! ASCII digits only (the Python's `isdigit` also took other Unicode
//! digits, which its `int()` then refused with a traceback).

use crate::y::{is_file, path_join};
use crate::ytree::{MIN_DEPTH, PAR1_END, PAR2_START, TABLES};
use crate::{Error, Outcome};
use og_filters::args::{self, Spec};
use og_io::Utf8;
use std::fmt::Write as _;
use std::io::Write;

const MIN_GQ: i128 = 20;
/// Yleaf's `-b` default (90%).
const MIN_FRACTION: f64 = 0.9;
const MAX_DEPTH_FACTOR: f64 = 2.0;
const MASK: &str = "callable_mask_chrY.fasta.gz";

/// The gates, in the order a record meets them and the `# dropped` line
/// prints them.
#[derive(Debug, Clone, Copy)]
enum Gate {
    Par,
    OutsideCallable,
    Filtered,
    Genotype,
    Depth,
    Quality,
    Fraction,
    ExcessDepth,
    Public,
}

const DROPPED: [(Gate, &str); 9] = [
    (Gate::Par, "PAR"),
    (Gate::OutsideCallable, "outside_callable"),
    (Gate::Filtered, "not_PASS"),
    (Gate::Genotype, "genotype"),
    (Gate::Depth, "depth"),
    (Gate::Quality, "quality"),
    (Gate::Fraction, "alt_fraction"),
    (Gate::ExcessDepth, "excess_depth"),
    (Gate::Public, "on_public_lists"),
];

/// The thresholds, as given.
struct Gates {
    min_depth: i128,
    min_gq: i128,
    min_fraction: f64,
    max_depth_factor: f64,
}

/// The chrY accessibility mask: one bit per base, set where the base is `P`.
#[derive(Debug)]
struct Mask {
    bits: Vec<u64>,
    len: u64,
}

impl Mask {
    fn passes(&self, pos: u64) -> bool {
        // 1-based; the Python indexed mask[pos - 1].
        pos >= 1
            && pos <= self.len
            && self.bits[((pos - 1) / 64) as usize] >> ((pos - 1) % 64) & 1 == 1
    }
}

fn read_mask(path: &str) -> Result<Mask, Error> {
    if !is_file(path) {
        return Err(Error(format!(
            "no chrY callable mask at {path}: install the current Y tree with scripts/setup.sh --yleaf-data <genome_dir>"
        )));
    }
    let mut m = Mask {
        bits: Vec::new(),
        len: 0,
    };
    let mut records = 0;
    let mut passed = 0u64;
    for line in crate::lines(path, Utf8::Strict)? {
        let line = line?;
        if line.starts_with('>') {
            records += 1;
            if records > 1 {
                return Err(Error(format!(
                    "the chrY callable mask {path} holds more than one sequence: it must be chrY alone"
                )));
            }
            continue;
        }
        // The Python joined the stripped lines and indexed the UTF-8 bytes.
        for &b in crate::strip(&line).as_bytes() {
            let bit = m.len % 64;
            if bit == 0 {
                m.bits.push(0);
            }
            if b == b'P'
                && let Some(w) = m.bits.last_mut()
            {
                *w |= 1 << bit;
                passed += 1;
            }
            m.len += 1;
        }
    }
    if passed == 0 {
        return Err(Error(format!(
            "the chrY callable mask {path} passes no base"
        )));
    }
    Ok(m)
}

/// Every position the two tables name, sorted and distinct.
fn public_positions(dir: &str) -> Result<Vec<u64>, Error> {
    let mut pos = Vec::new();
    let mut found = 0;
    for t in TABLES {
        let path = path_join(dir, t);
        if !is_file(&path) {
            continue;
        }
        found += 1;
        for line in crate::lines(&path, Utf8::Strict)? {
            let line = line?;
            let mut c = line.strip_suffix('\n').unwrap_or(&line).split('\t');
            if let Some(p) = c.nth(3).filter(|p| crate::digits(p)) {
                // A position past u64 is past the mask's end: no record reaches it.
                if let Ok(p) = p.parse::<u64>() {
                    pos.push(p);
                }
            }
        }
    }
    if found < TABLES.len() {
        return Err(Error(format!(
            "{dir} has no {}: install the current Y tree with scripts/setup.sh --yleaf-data <genome_dir>",
            TABLES.join(" and no ")
        )));
    }
    if pos.is_empty() {
        return Err(Error(format!("the SNP tables in {dir} name no position")));
    }
    pos.sort_unstable();
    pos.dedup();
    Ok(pos)
}

/// `FORMAT:sample` lookup: the value of the first `key`, `None` when the key
/// or its value is missing.
fn field<'a>(fmt: &str, sample: &'a str, key: &str) -> Option<&'a str> {
    let i = fmt.split(':').position(|k| k == key)?;
    sample.split(':').nth(i)
}

/// The genotype is the ALT allele alone.
fn carries_alt(gt: Option<&str>) -> bool {
    gt.is_some_and(|gt| gt.split(['/', '|']).all(|a| a == "1"))
}

/// `int(s)` for ASCII digits.
fn uint(s: &str) -> Option<u128> {
    if crate::digits(s) {
        s.parse().ok()
    } else {
        None
    }
}

/// (depth, ALT fraction) from AD when it is usable, else from DP and VAF.
fn depth_and_fraction(fmt: &str, sample: &str) -> (Option<u128>, Option<f64>) {
    if let Some(ad) = field(fmt, sample, "AD").filter(|s| !s.is_empty()) {
        let counts: Option<Vec<u128>> = ad.split(',').map(uint).collect();
        if let Some(counts) = counts.filter(|c| c.len() >= 2) {
            let total: u128 = counts.iter().sum();
            if total > 0 {
                return (Some(total), Some(counts[1] as f64 / total as f64));
            }
        }
    }
    match field(fmt, sample, "DP").and_then(uint) {
        Some(dp) if dp > 0 => (Some(dp), field(fmt, sample, "VAF").and_then(crate::float)),
        _ => (None, None),
    }
}

/// A call that passed gates 1 to 7.
struct Kept {
    pos: u64,
    ref_: char,
    alt: char,
    depth: u128,
    fraction: Option<f64>,
    gq: u128,
    gt: String,
}

/// The median the Python's `median()` gave: an int for an odd count, the
/// float mean of the middle two for an even one.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Median {
    Int(u128),
    Float(f64),
}

impl Median {
    fn of(depths: &mut [u128]) -> Option<Median> {
        depths.sort_unstable();
        let n = depths.len();
        match n {
            0 => None,
            _ if n % 2 == 1 => Some(Median::Int(depths[n / 2])),
            _ => Some(Median::Float(
                (depths[n / 2 - 1] + depths[n / 2]) as f64 / 2.0,
            )),
        }
    }

    fn value(self) -> f64 {
        match self {
            Median::Int(v) => v as f64,
            Median::Float(v) => v,
        }
    }

    /// `str()` of it.
    fn text(self) -> String {
        match self {
            Median::Int(v) => v.to_string(),
            Median::Float(v) => pycompat::fmt::repr(v),
        }
    }
}

/// One chrY SNV record before the gates, as far as the gates read it.
fn gate_record(c: &[&str], mask: &Mask, g: &Gates, drops: &mut [u64; 9]) -> Option<Kept> {
    let pos: u64 = c[1].parse().unwrap_or(u64::MAX);
    let fmt = c[8];
    let sample = c[9];
    let drop = |drops: &mut [u64; 9], gate: Gate| {
        drops[gate as usize] += 1;
        None
    };
    if pos <= PAR1_END || pos >= PAR2_START {
        return drop(drops, Gate::Par);
    }
    if !mask.passes(pos) {
        return drop(drops, Gate::OutsideCallable);
    }
    if !matches!(c[6], "PASS" | "." | "") {
        return drop(drops, Gate::Filtered);
    }
    let gt = field(fmt, sample, "GT");
    if !carries_alt(gt) {
        return drop(drops, Gate::Genotype);
    }
    let (depth, fraction) = depth_and_fraction(fmt, sample);
    let Some(depth) = depth.filter(|&d| d as i128 >= g.min_depth || d > i128::MAX as u128) else {
        return drop(drops, Gate::Depth);
    };
    let Some(gq) = field(fmt, sample, "GQ")
        .and_then(uint)
        .filter(|&q| q > i128::MAX as u128 || q as i128 >= g.min_gq)
    else {
        return drop(drops, Gate::Quality);
    };
    // NaN passes, as `nan < x` is false in Python.
    if fraction.is_none_or(|f| f < g.min_fraction) {
        return drop(drops, Gate::Fraction);
    }
    let base = |s: &str| {
        let mut ch = s.chars();
        ch.next().map(|c| c.to_ascii_uppercase())
    };
    Some(Kept {
        pos,
        ref_: base(c[3]).unwrap_or('N'),
        alt: base(c[4]).unwrap_or('N'),
        depth,
        fraction,
        gq,
        gt: gt.unwrap_or_default().to_owned(),
    })
}

/// One REF base and one ALT base, upper-cased as Python's `str.upper()`
/// would make them one of A, C, G, T.
fn is_snv_base(s: &str) -> bool {
    matches!(s, "A" | "C" | "G" | "T" | "a" | "c" | "g" | "t")
}

struct Run<'a> {
    vcf: &'a str,
    dir: &'a str,
    sample: &'a str,
    prefix: &'a str,
    mask: Option<&'a str>,
    gates: Gates,
}

fn private_snps(r: &Run<'_>) -> Outcome {
    let fail = |e| (e, 1);
    let public = public_positions(r.dir).map_err(fail)?;
    let mask_path = match r.mask {
        Some(m) => m.to_owned(),
        None => os_path_join(r.dir, MASK),
    };
    let mask = read_mask(&mask_path).map_err(fail)?;
    let mut drops = [0u64; 9];
    let mut records = 0u64;
    let mut kept: Vec<Kept> = Vec::new();
    for line in crate::lines(r.vcf, Utf8::Strict).map_err(fail)? {
        let line = line.map_err(fail)?;
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(&line)
            .split('\t')
            .collect();
        if c.len() < 10 || !matches!(c[0], "chrY" | "Y") || !crate::digits(c[1]) {
            continue;
        }
        if !is_snv_base(c[3]) || !is_snv_base(c[4]) {
            continue;
        }
        records += 1;
        if let Some(k) = gate_record(&c, &mask, &r.gates, &mut drops) {
            kept.push(k);
        }
    }
    // The single-copy gate needs the sample's own chrY depth: the median of
    // the calls that passed the rest.
    let median = Median::of(&mut kept.iter().map(|k| k.depth).collect::<Vec<_>>());
    let limit = median
        .filter(|m| m.value() != 0.0 && r.gates.max_depth_factor > 0.0)
        .map(|m| m.value() * r.gates.max_depth_factor);
    let mut rows: Vec<Kept> = Vec::new();
    for k in kept {
        if limit.is_some_and(|l| k.depth as f64 > l) {
            drops[Gate::ExcessDepth as usize] += 1;
        } else if public.binary_search(&k.pos).is_ok() {
            drops[Gate::Public as usize] += 1;
        } else {
            rows.push(k);
        }
    }
    rows.sort_by_key(|k| k.pos);
    let median_text = median.map_or_else(|| "0".to_owned(), Median::text);
    let mask_name = basename(r.mask.unwrap_or(MASK));
    let g = &r.gates;
    let mut body = format!(
        "# sample\t{}\n# gates\tmin_depth={} min_gq={} min_fraction={} max_depth_factor={} \
         median_chrY_depth={median_text} callable_mask={mask_name}\n# dropped\t",
        r.sample,
        g.min_depth,
        g.min_gq,
        pycompat::fmt::repr(g.min_fraction),
        pycompat::fmt::repr(g.max_depth_factor),
    );
    let dropped: Vec<String> = DROPPED
        .iter()
        .map(|(gate, name)| format!("{name}={}", drops[*gate as usize]))
        .collect();
    body.push_str(&dropped.join(" "));
    body.push_str("\nPos\tRef\tAlt\tDepth\tAlt_fraction\tGQ\tGenotype\n");
    for k in &rows {
        let frac = k
            .fraction
            .map_or_else(|| ".".to_owned(), |f| crate::fixed(f, 3));
        let _ = writeln!(
            body,
            "{}\t{}\t{}\t{}\t{frac}\t{}\t{}",
            k.pos, k.ref_, k.alt, k.depth, k.gq, k.gt
        );
    }
    crate::write_atomic(&format!("{}_y_private_snps.tsv", r.prefix), body.as_bytes())
        .map_err(fail)?;
    Ok(format!(
        "Private chrY SNVs: {} past the public lists ({records} chrY SNV records read, {} outside the callable mask, \
         {} already named by a public list, median chrY depth {median_text})\n  Each one is unconfirmed until \
         another man is found who carries it; some are read errors (docs/37-y-haplogroup.md).\n",
        rows.len(),
        drops[Gate::OutsideCallable as usize],
        drops[Gate::Public as usize],
    ))
}

/// `os.path.join(dir, name)` for a relative `name`.
fn os_path_join(dir: &str, name: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') {
        format!("{dir}{name}")
    } else {
        format!("{dir}/{name}")
    }
}

/// `os.path.basename`.
fn basename(p: &str) -> &str {
    p.rsplit('/').next().unwrap_or(p)
}

const SPEC: Spec = Spec {
    prog: "og y private-snps",
    usage: "usage: og y private-snps [-h] --vcf VCF --tree-dir TREE_DIR --sample SAMPLE --out-prefix OUT_PREFIX \
[--callable-mask CALLABLE_MASK] [--min-depth MIN_DEPTH] [--min-gq MIN_GQ] [--min-fraction MIN_FRACTION] \
[--max-depth-factor MAX_DEPTH_FACTOR]\n",
    options: &[
        "--vcf",
        "--tree-dir",
        "--sample",
        "--out-prefix",
        "--callable-mask",
        "--min-depth",
        "--min-gq",
        "--min-fraction",
        "--max-depth-factor",
    ],
    positionals: &[],
    required: &["--vcf", "--tree-dir", "--sample", "--out-prefix"],
};

pub(crate) fn cli(rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    crate::run(&SPEC, rest, out, err, |a| {
        let gates = Gates {
            min_depth: a.typed("--min-depth", i128::from(MIN_DEPTH), "int", args::int)?,
            min_gq: a.typed("--min-gq", MIN_GQ, "int", args::int)?,
            min_fraction: a.typed("--min-fraction", MIN_FRACTION, "float", args::float)?,
            max_depth_factor: a.typed(
                "--max-depth-factor",
                MAX_DEPTH_FACTOR,
                "float",
                args::float,
            )?,
        };
        let refuse = |m: &str| Ok(Err((Error(m.to_owned()), 1)));
        if gates.min_depth < 1 {
            return refuse("--min-depth must be 1 or more");
        }
        // `not 0 < f <= 1`, so NaN is refused.
        if !(gates.min_fraction > 0.0 && gates.min_fraction <= 1.0) {
            return refuse("--min-fraction must be above 0 and at most 1");
        }
        // NaN is refused too: the Python's `< 0` let it through and turned
        // the depth gate off.
        if gates.max_depth_factor.is_nan() || gates.max_depth_factor < 0.0 {
            return refuse("--max-depth-factor must be 0 or more");
        }
        Ok(private_snps(&Run {
            vcf: crate::opt(a, "--vcf"),
            dir: crate::opt(a, "--tree-dir"),
            sample: crate::opt(a, "--sample"),
            prefix: crate::opt(a, "--out-prefix"),
            // An empty --callable-mask is the default, as Python's `or`.
            mask: a.get("--callable-mask").filter(|m| !m.is_empty()),
            gates,
        }))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    #[test]
    fn median_keeps_the_int_or_float() {
        assert_eq!(
            Median::of(&mut [31, 5, 40]).map(Median::text),
            Some("31".into())
        );
        assert_eq!(
            Median::of(&mut [31, 32]).map(Median::text),
            Some("31.5".into())
        );
        assert_eq!(
            Median::of(&mut [30, 32]).map(Median::text),
            Some("31.0".into())
        );
        assert_eq!(Median::of(&mut []), None);
    }

    #[test]
    fn format_fields() {
        assert_eq!(field("GT:AD:DP", "1:0,30:30", "AD"), Some("0,30"));
        assert_eq!(field("GT:AD:DP", "1:0,30", "DP"), None);
        assert_eq!(field("GT:AD", "1:0,30", "GQ"), None);
        assert!(carries_alt(Some("1|1")) && carries_alt(Some("1")));
        assert!(!carries_alt(Some("0/1")) && !carries_alt(Some(".")) && !carries_alt(None));
        assert_eq!(depth_and_fraction("AD", "6,24"), (Some(30), Some(0.8)));
        assert_eq!(
            depth_and_fraction("AD:DP:VAF", "0,0:12:0.95"),
            (Some(12), Some(0.95))
        );
        assert_eq!(depth_and_fraction("AD:DP:VAF", "x:12:no"), (Some(12), None));
        assert_eq!(depth_and_fraction("AD:DP", "5:0"), (None, None));
    }

    fn mask_file(d: &Dir, name: &str, text: &str) -> String {
        d.put(name, text)
    }

    /// Red first: indexing the mask at `pos` instead of `pos - 1` would pass
    /// 3 and drop 1.
    #[test]
    fn red_first_mask_is_one_based() {
        let d = Dir::new("yprivate-mask");
        let m = read_mask(&mask_file(&d, "m.fa", ">chrY\nPZ\nP P\n")).unwrap();
        assert_eq!(m.len, 5, "stripped lines are joined");
        assert_eq!(
            (1..=6).map(|p| m.passes(p)).collect::<Vec<_>>(),
            [true, false, true, false, true, false]
        );
        assert!(!m.passes(0));
        let two = mask_file(&d, "two.fa", ">chrY\nPP\n>chrX\nPP\n");
        assert!(
            read_mask(&two)
                .unwrap_err()
                .0
                .contains("more than one sequence")
        );
        let none = mask_file(&d, "none.fa", ">chrY\nNN\n");
        assert!(read_mask(&none).unwrap_err().0.contains("passes no base"));
    }

    #[test]
    fn a_long_mask_spans_words() {
        let d = Dir::new("yprivate-long");
        let seq: String = (0..200)
            .map(|i| if i % 3 == 0 { 'P' } else { 'N' })
            .collect();
        let m = read_mask(&mask_file(&d, "m.fa", &format!(">chrY\n{seq}\n"))).unwrap();
        assert!((1..=200u64).all(|p| m.passes(p) == ((p - 1) % 3 == 0)));
    }

    /// A NaN depth factor would turn the depth gate off without a word.
    #[test]
    fn a_nan_depth_factor_is_refused() {
        for (factor, refused) in [("nan", true), ("-1", true), ("0", false)] {
            let args: Vec<String> = [
                "--vcf",
                "/nonexistent.vcf",
                "--tree-dir",
                "/nonexistent",
                "--sample",
                "S",
                "--out-prefix",
                "/nonexistent/o",
                "--max-depth-factor",
                factor,
            ]
            .map(String::from)
            .into();
            let (mut out, mut err) = (Vec::new(), Vec::new());
            assert_eq!(cli(&args, &mut out, &mut err), 1, "{factor}");
            let err = String::from_utf8(err).unwrap();
            assert_eq!(
                err.contains("--max-depth-factor must be 0 or more"),
                refused,
                "{factor}: {err}"
            );
        }
    }
}
