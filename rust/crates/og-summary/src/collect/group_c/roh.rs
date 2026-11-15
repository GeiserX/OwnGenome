//! Step 11 (bcftools roh): the runs of homozygosity of chr1 to chr22
//! outside GRCh38's centromere and heterochromatin bands, F_ROH, the
//! relationship band it is nearest to and the one-chromosome (UPD) flag;
//! and `og roh summary`, the step's own table of the long runs.
//!
//! [`MASKED_BANDS`] is the one copy of the bands: the summary section,
//! CNVpytor's count and `og roh summary` (step 11 and the ROH module) all
//! read it here.

use super::util::{Obj, count, lines, s};
use super::{Ctx, Out};
use pycompat::json::Value;
use std::io::Write;

/// GRCh38's centromeres and heterochromatin: the acen, gvar and stalk rows
/// of UCSC's cytoBand.hg38.txt (the file setup.sh installs as
/// reference/cytoBand.hg38.txt), neighbouring rows merged, as (chrom,
/// 0-based start, exclusive end). Reads there map badly, so a run of
/// homozygosity or a copy number call over them is mostly an artefact.
/// Regenerate with
///   gzip -cd cytoBand.txt.gz | awk '$1 ~ /^chr([0-9]+|X|Y)$/ && $5 ~ /^(acen|gvar|stalk)$/' \
///     | sort -k1,1V -k2,2n | awk '{ if ($1 == c && $2 == e) e = $3; else { if (c != "") printf "%s:%d-%d ", c, s, e;
///       c = $1; s = $2; e = $3 } } END { print c ":" s "-" e }'
pub const MASKED_BANDS: [(&str, i128, i128); 25] = [
    ("chr1", 121_700_000, 143_200_000),
    ("chr2", 91_800_000, 96_000_000),
    ("chr3", 87_800_000, 98_600_000),
    ("chr4", 48_200_000, 51_800_000),
    ("chr5", 46_100_000, 51_400_000),
    ("chr6", 58_500_000, 62_600_000),
    ("chr7", 58_100_000, 62_100_000),
    ("chr8", 43_200_000, 47_200_000),
    ("chr9", 42_200_000, 61_500_000),
    ("chr10", 38_000_000, 41_600_000),
    ("chr11", 51_000_000, 55_800_000),
    ("chr12", 33_200_000, 37_800_000),
    ("chr13", 0, 18_900_000),
    ("chr14", 0, 18_200_000),
    ("chr15", 0, 20_500_000),
    ("chr16", 35_300_000, 47_000_000),
    ("chr17", 22_700_000, 27_400_000),
    ("chr18", 15_400_000, 21_500_000),
    ("chr19", 19_900_000, 31_900_000),
    ("chr20", 25_700_000, 30_400_000),
    ("chr21", 0, 13_000_000),
    ("chr22", 0, 17_400_000),
    ("chrX", 58_100_000, 63_800_000),
    ("chrY", 10_300_000, 10_600_000),
    ("chrY", 26_600_000, 57_227_415),
];

/// GRCh38's chr1 to chr22 (UCSC hg38.chrom.sizes).
pub const AUTOSOME_BP: i128 = 2_875_001_522;

// UCSC hg38.chrom.sizes, in chromosome order. ST rows tell us which
// chromosome lengths belong in this run's denominator, not callable bases.
const AUTOSOME_LENGTHS: [i128; 22] = [
    248_956_422,
    242_193_529,
    198_295_559,
    190_214_555,
    181_538_259,
    170_805_979,
    159_345_973,
    145_138_636,
    138_394_717,
    133_797_422,
    135_086_622,
    133_275_309,
    114_364_328,
    107_043_718,
    101_991_189,
    90_338_345,
    83_257_441,
    80_373_285,
    58_617_616,
    64_444_167,
    46_709_983,
    50_818_468,
];

fn autosome_index(chrom: &str) -> Option<usize> {
    is_autosome(chrom).then(|| chrom[3..].parse::<usize>().unwrap() - 1)
}

fn covered_unmasked_bp(covered: &[bool; 22]) -> i128 {
    covered
        .iter()
        .enumerate()
        .filter(|(_, yes)| **yes)
        .map(|(i, _)| {
            let chrom = format!("chr{}", i + 1);
            AUTOSOME_LENGTHS[i]
                - MASKED_BANDS
                    .iter()
                    .filter(|(c, ..)| *c == chrom)
                    .map(|(_, a, b)| b - a)
                    .sum::<i128>()
        })
        .sum()
}
/// The runs F_ROH and the one-chromosome flag count (docs/11-roh-analysis.md).
pub const ROH_MIN_BP: f64 = 5e6;
/// One chromosome carrying most of the long runs: at least this much of it
/// in runs of 5 Mb or more ...
pub const UPD_MIN_BP: f64 = 20e6;
/// ... and more than this share of all such autosomal runs.
pub const UPD_SHARE: f64 = 0.5;

/// Bases of the 1-based, inclusive span `start..=end` inside the masked
/// bands of `chrom`.
pub fn masked_overlap(chrom: &str, start: i128, end: i128) -> i128 {
    MASKED_BANDS
        .iter()
        .filter(|(c, ..)| *c == chrom)
        .map(|(_, a, b)| overlap(*a, *b, start, end))
        .sum()
}

fn overlap(a: i128, b: i128, start: i128, end: i128) -> i128 {
    (end.min(b) - (start - 1).max(a)).max(0)
}

/// Whether a chromosome name is one of chr1 to chr22.
pub fn is_autosome(chrom: &str) -> bool {
    chrom
        .strip_prefix("chr")
        .and_then(|n| {
            n.parse::<u8>()
                .ok()
                .filter(|v| (1..=22).contains(v))
                .map(|v| (v, n))
        })
        .is_some_and(|(v, n)| n == v.to_string())
}

/// [`AUTOSOME_BP`] without the autosomal bands: the length F_ROH is a share
/// of, since a run's part in a band is not counted either.
pub fn autosome_unmasked_bp() -> i128 {
    AUTOSOME_BP
        - MASKED_BANDS
            .iter()
            .filter(|(c, ..)| is_autosome(c))
            .map(|(_, a, b)| b - a)
            .sum::<i128>()
}

/// docs/11-roh-analysis.md's table: F of a child of the parents. A value
/// goes to the row whose F is nearest on a log scale (the boundaries are
/// the geometric means of neighbouring rows), and under 1/128 (half the
/// second-cousin F) to "no sign of related parents".
pub fn roh_band(f: f64) -> &'static str {
    let gm = |a: f64, b: f64| (a * b).powf(0.5);
    let bands = [
        (1.0 / 128.0, "no sign of related parents"),
        (gm(1.0 / 64.0, 1.0 / 32.0), "second cousins (F 1/64)"),
        (
            gm(1.0 / 32.0, 1.0 / 16.0),
            "first cousins once removed (F 1/32)",
        ),
        (gm(1.0 / 16.0, 1.0 / 8.0), "first cousins (F 1/16)"),
        (
            gm(1.0 / 8.0, 1.0 / 4.0),
            "half siblings, uncle and niece, double first cousins (F 1/8)",
        ),
        (f64::INFINITY, "parent and child, full siblings (F 1/4)"),
    ];
    bands
        .iter()
        .find(|(top, _)| f < *top)
        .map_or("parent and child, full siblings (F 1/4)", |b| b.1)
}

/// One source for region names and their citations. Embed the shipped
/// asset so a staged report needs no second runtime file or copied table.
const IMPRINTED_TSV: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../assets/imprinted_regions.tsv"
));

/// A chromosome hint cannot determine parental origin or diagnose a
/// syndrome. The two legacy keys stay as placeholders for old readers.
fn imprinted(chrom: &str) -> Value {
    Value::Array(
        IMPRINTED_TSV
            .lines()
            .filter(|line| !line.starts_with('#'))
            .skip(1)
            .filter_map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                (f.len() == 7 && f[0] == chrom).then(|| {
                    Obj::new()
                        .set("region", s(f[3]))
                        .set("maternal_upd", s("-"))
                        .set("paternal_upd", s("-"))
                        .set("pmid", s(f[6]))
                        .val()
                })
            })
            .collect(),
    )
}

fn r1(x: f64, n: i32) -> Value {
    super::util::round(x, n)
}

/// Step 11's bcftools roh output: runs on chrX and chrY counted apart, runs
/// on any contig other than chr1 to chr22 left out, and the part of a run
/// in a masked band not counted.
pub fn section(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find(&[
        format!("vcf/{s_}_roh.txt"),
        format!("roh/{s_}_roh.txt"),
        "roh/*_roh.txt".to_owned(),
    ]) else {
        return Ok(None);
    };
    let (mut total, mut largest, mut masked) = (0.0_f64, 0.0_f64, 0.0_f64);
    let (mut segs, mut sex_segs) = (0usize, 0usize);
    let mut big = Vec::new();
    // chrom -> the bases of its long runs
    let mut per_chrom: indexmap::IndexMap<String, f64> = indexmap::IndexMap::new();
    let (mut assessed, mut observed) = ([false; 22], [false; 22]);
    let mut reader = crate::text_lines(&p.path)?;
    let mut line = String::new();
    loop {
        line.clear();
        if !reader.read_line(&mut line).map_err(|e| e.to_string())? {
            break;
        }
        let c: Vec<&str> = line.split_whitespace().collect();
        if c.len() >= 6
            && c[0] == "ST"
            && let Some(i) = autosome_index(c[2])
            && crate::int(c[3]).is_ok_and(|pos| pos > 0 && pos <= AUTOSOME_LENGTHS[i])
            && matches!(c[4], "0" | "1")
        {
            assessed[i] = true;
        }
        if c.len() < 6 || c[0] != "RG" {
            continue;
        }
        let (Ok(start), Ok(end), Ok(mut length)) = (
            crate::int(c[3]),
            crate::int(c[4]),
            pycompat::num::parse_float(c[5]),
        ) else {
            continue;
        };
        if !length.is_finite() {
            return Err(format!(
                "{} has a run of length {}",
                crate::basename(&p.path),
                c[5]
            ));
        }
        if matches!(c[2], "chrX" | "chrY" | "X" | "Y") {
            sex_segs += 1;
            continue;
        }
        if !is_autosome(c[2]) {
            continue;
        }
        observed[autosome_index(c[2]).unwrap()] = true;
        let inside = length.min(masked_overlap(c[2], start, end) as f64);
        masked += inside;
        length -= inside;
        if length <= 0.0 {
            continue;
        }
        segs += 1;
        total += length;
        largest = largest.max(length);
        if length >= ROH_MIN_BP {
            let pc = per_chrom.entry(c[2].to_owned()).or_default();
            *pc += length;
            let mut row = Obj::new()
                .set("region", s(&format!("{}:{}-{}", c[2], c[3], c[4])))
                .set("mb", r1(length / 1e6, 1));
            if inside != 0.0 {
                row.put("masked_mb", r1(inside / 1e6, 1));
            }
            big.push(row.val());
        }
    }
    let long_bp: f64 = per_chrom.values().sum();
    let has_sites = assessed.iter().any(|yes| *yes);
    // Old -O r files omit ST rows. Their runs give a lower bound on
    // chromosome scope, never evidence for a genome-wide interpretation.
    let covered: [bool; 22] = std::array::from_fn(|i| assessed[i] || observed[i]);
    let denominator = covered_unmasked_bp(&covered);
    let complete = has_sites
        && assessed.iter().all(|yes| *yes)
        && denominator >= autosome_unmasked_bp()
        && observed
            .iter()
            .zip(assessed)
            .all(|(run, site)| !run || site);
    let f_roh = (denominator > 0).then(|| long_bp / denominator as f64);
    let mut upd = Value::Null;
    // Retain the first chromosome on ties, as before.
    let top = per_chrom
        .iter()
        .fold(None::<(&String, &f64)>, |best, kv| match best {
            Some(b) if b.1 >= kv.1 => Some(b),
            _ => Some(kv),
        });
    if complete
        && let Some((chrom, bp)) = top
        && *bp >= UPD_MIN_BP
        && *bp > UPD_SHARE * long_bp
    {
        upd = Obj::new()
            .set("chrom", s(chrom))
            .set("mb", r1(bp / 1e6, 1))
            .set("share", r1(bp / long_bp, 2))
            .set("imprinted", imprinted(chrom))
            .set("note", s("Long chromosome-focused ROH; investigate possible isodisomy. A deletion, acquired loss of heterozygosity or shared ancestry can give this pattern. Copy neutrality and parental origin are not established; check copy-number results and parents' genotypes. Heterodisomy can leave no ROH."))
            .val();
    }
    let relationship = if complete && upd == Value::Null {
        s(roh_band(f_roh.unwrap_or(0.0)))
    } else {
        Value::Null
    };
    let data = Obj::new()
        .set("segments", count(segs))
        .set("total_mb", r1(total / 1e6, 1))
        .set("largest_mb", r1(largest / 1e6, 1))
        .set("autosomal_over_5mb", Value::Array(big))
        .set("masked_mb", r1(masked / 1e6, 1))
        .set("sex_chromosome_segments", count(sex_segs))
        .set("long_mb", r1(long_bp / 1e6, 1))
        .set("f_roh", f_roh.map_or(Value::Null, |f| r1(f, 4)))
        .set("relationship", relationship)
        .set("upd_flag", upd)
        .set("scope", Obj::new()
            .set("source", s(if has_sites { "sites" } else { "segments_only" }))
            .set("autosomes", Value::Array(covered.iter().enumerate()
                .filter(|(_, yes)| **yes).map(|(i, _)| s(&format!("chr{}", i + 1))).collect()))
            .set("unmasked_bp", Value::Int(denominator))
            .set("partial", Value::Bool(!complete))
            .set("note", s(if complete {
                "Sites on all 22 autosomes; marker presence does not establish callable breadth or copy neutrality."
            } else {
                "Partial or unverified chromosome scope: no relationship band or UPD flag. F_ROH uses only the listed chromosomes, not callable bases."
            })).val());
    Ok(Some((p, data.val())))
}

// --- og roh summary ------------------------------------------------------------

const USAGE: &str = "usage: og roh summary --sample SAMPLE --roh ROH [--out OUT]\n";

/// A field as awk reads it as a number: the longest leading decimal
/// number, else 0.
fn awk_num(x: &str) -> f64 {
    let b = x.as_bytes();
    let mut i = 0;
    if matches!(b.first(), Some(b'+' | b'-')) {
        i += 1;
    }
    let digits = |i: &mut usize| {
        let s = *i;
        while b.get(*i).is_some_and(u8::is_ascii_digit) {
            *i += 1;
        }
        *i - s
    };
    let mut n = digits(&mut i);
    if b.get(i) == Some(&b'.') {
        i += 1;
        n += digits(&mut i);
    }
    if n == 0 {
        return 0.0;
    }
    let mantissa_end = i;
    if matches!(b.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(b.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        if digits(&mut j) > 0 {
            i = j;
        } else {
            i = mantissa_end;
        }
    }
    x[..i].parse().unwrap_or(0.0)
}

/// The rows of step 11's summary: every run on chr1 to chr22 whose part
/// outside the masked bands is 5 Mb or more, with that part's length.
pub fn summary_rows(roh: &str) -> Result<String, String> {
    let mut out = String::new();
    for line in lines(roh)? {
        // awk's fields: runs of blanks, leading ones dropped.
        let f: Vec<&str> = line.split([' ', '\t']).filter(|x| !x.is_empty()).collect();
        if f.len() < 3 || f[0] != "RG" || !is_autosome(f[2]) {
            continue;
        }
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let (start, end, len) = (awk_num(get(3)), awk_num(get(4)), awk_num(get(5)));
        let mut m = 0.0;
        for (c, a, b) in MASKED_BANDS {
            if c == f[2] {
                let o = end.min(b as f64) - (start - 1.0).max(a as f64);
                if o > 0.0 {
                    m += o;
                }
            }
        }
        if m > len {
            m = len;
        }
        let l = len - m;
        if l >= 5_000_000.0 {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                f[2],
                get(3),
                get(4),
                l.trunc() as i64,
                pycompat::fmt::fixed(l / 1e6, 1)
            ));
        }
    }
    Ok(out)
}

/// `og roh summary`: step 11's summary table of the long runs, for the bash
/// step and the ROH module alike, into `--out` or on stdout. Exit 0 when
/// written, 1 when the roh file cannot be read or the table written, 2 for
/// a command line it cannot read.
pub fn cli(argv: &[String], out: &mut impl Write, err: &mut impl Write) -> u8 {
    use crate::args::{Kind, Opt, parse};
    const OPTS: &[Opt] = &[
        Opt("--sample", Kind::Value, true),
        Opt("--roh", Kind::Value, true),
        Opt("--out", Kind::Value, false),
    ];
    let a = match parse(argv, OPTS) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og roh summary: error: {e}");
            return 2;
        }
    };
    let mut run = || -> Result<(), String> {
        let rows = summary_rows(a.text("--roh"))?;
        let body = format!(
            "# ROH Summary for {}\n# Segments >=5MB on autosomes (excludes chrX/chrY)\nchrom\tstart\tend\tlength_bp\tlength_mb\n{rows}",
            a.text("--sample")
        );
        match a.get("--out") {
            Some(path) => crate::write_file(path, body.as_bytes()),
            None => out
                .write_all(body.as_bytes())
                .and_then(|()| out.flush())
                .map_err(|e| format!("stdout: {e}")),
        }
    };
    crate::finish(run(), err)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    #[test]
    fn imprinted_asset_is_the_only_table() {
        let rows: Vec<_> = IMPRINTED_TSV
            .lines()
            .filter(|l| !l.starts_with('#'))
            .skip(1)
            .collect();
        assert_eq!(rows.len(), 6);
        for chrom in ["chr6", "chr7", "chr11", "chr14", "chr15", "chr20"] {
            let Value::Array(hints) = imprinted(chrom) else {
                panic!()
            };
            assert_eq!(hints.len(), 1);
            let row = rows
                .iter()
                .find(|line| line.starts_with(&format!("{chrom}\t")))
                .unwrap();
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(f.len(), 7);
            assert!(f[6].split(';').all(|pmid| pmid.parse::<u32>().is_ok()));
            assert!(f[1].is_empty() && f[2].is_empty());
        }
        assert_eq!(AUTOSOME_LENGTHS.iter().sum::<i128>(), AUTOSOME_BP);
        assert_eq!(covered_unmasked_bp(&[true; 22]), autosome_unmasked_bp());
    }

    #[test]
    fn awk_reads_the_leading_number() {
        assert_eq!(awk_num("7000000"), 7_000_000.0);
        assert_eq!(awk_num("8000000.5"), 8_000_000.5);
        assert_eq!(awk_num("1e3x"), 1000.0);
        assert_eq!(awk_num("1e"), 1.0);
        assert_eq!(awk_num("x"), 0.0);
        assert_eq!(awk_num("-2"), -2.0);
    }

    #[test]
    fn the_summary_keeps_the_unmasked_part_of_long_autosomal_runs() {
        let d = Dir::new("roh-summary");
        let roh = d.put(
            "roh.txt",
            "# RG\nRG\tS\tchrY\t1\t9000000\t9000000\t9\t9\n\
             RG\tS\tchr1\t120000001\t145000000\t25000000\t9\t9\n\
             RG\tS\tchr2\t10000001\t17000000\t7000000\t9\t9\n\
             RG\tS\tchr13\t1\t25000000\t25000000\t9\t9\n\
             RG\tS\tchr5\t1\t4999999\t4999999\t9\t9\n",
        );
        let argv: Vec<String> = ["--sample", "S", "--roh", &roh]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(cli(&argv, &mut out, &mut err), 0);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "# ROH Summary for S\n# Segments >=5MB on autosomes (excludes chrX/chrY)\n\
             chrom\tstart\tend\tlength_bp\tlength_mb\n\
             chr2\t10000001\t17000000\t7000000\t7.0\n\
             chr13\t1\t25000000\t6100000\t6.1\n"
        );
    }

    #[test]
    fn a_missing_roh_file_exits_1_and_no_sample_exits_2() {
        let argv = |a: &[&str]| a.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(
            cli(
                &argv(&["--sample", "S", "--roh", "/no/such"]),
                &mut out,
                &mut err
            ),
            1
        );
        assert_eq!(cli(&argv(&["--roh", "x"]), &mut out, &mut err), 2);
    }
}
