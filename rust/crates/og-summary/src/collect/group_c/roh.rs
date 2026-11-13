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

/// assets/imprinted_regions.tsv, the same rows (a test keeps the two
/// equal): chrom, the region's span (none for the whole chromosome), the
/// region, the maternal and paternal UPD syndromes and the PMIDs.
/// chrom, span, then four text fields.
type Imprinted = (
    &'static str,
    Option<(i128, i128)>,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
);

pub const IMPRINTED_REGIONS: [Imprinted; 6] = [
    (
        "chr6",
        None,
        "6q24 (PLAGL1, HYMAI)",
        "-",
        "Transient neonatal diabetes mellitus (6q24-related)",
        "20301706",
    ),
    (
        "chr7",
        None,
        "7p12 (GRB10), 7q32 (MEST)",
        "Silver-Russell syndrome",
        "-",
        "27585961",
    ),
    (
        "chr11",
        Some((0, 22_000_000)),
        "11p15 (H19/IGF2, KCNQ1OT1)",
        "Silver-Russell syndrome",
        "Beckwith-Wiedemann syndrome",
        "29377879;27585961",
    ),
    (
        "chr14",
        None,
        "14q32.2 (DLK1-MEG3)",
        "Temple syndrome",
        "Kagami-Ogata syndrome",
        "24891339;26377239",
    ),
    (
        "chr15",
        Some((20_500_000, 33_400_000)),
        "15q11-q13 (SNRPN, UBE3A)",
        "Prader-Willi syndrome",
        "Angelman syndrome",
        "22237428;27615419",
    ),
    (
        "chr20",
        None,
        "20q13.32 (GNAS)",
        "Mulchandani-Bhoj-Conlin syndrome (growth failure)",
        "Pseudohypoparathyroidism type 1B",
        "26248010;29959430",
    ),
];

/// The imprinted regions a flagged chromosome joins: every row of the
/// whole chromosome, and a row with a region only when a run overlaps it.
fn imprinted(chrom: &str, runs: &[(i128, i128)]) -> Value {
    Value::Array(
        IMPRINTED_REGIONS
            .iter()
            .filter(|r| r.0 == chrom)
            .filter(|r| {
                r.1.is_none_or(|(a, b)| runs.iter().any(|&(s, e)| overlap(a, b, s, e) > 0))
            })
            .map(|r| {
                Obj::new()
                    .set("region", s(r.2))
                    .set("maternal_upd", s(r.3))
                    .set("paternal_upd", s(r.4))
                    .set("pmid", s(r.5))
                    .val()
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
    // chrom -> the bases of its long runs, and the runs (start, end)
    let mut per_chrom: indexmap::IndexMap<String, (f64, Vec<(i128, i128)>)> =
        indexmap::IndexMap::new();
    for line in lines(&p.path)? {
        let c = pycompat::text::split(&line, None, None);
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
            pc.0 += length;
            pc.1.push((start, end));
            let mut row = Obj::new()
                .set("region", s(&format!("{}:{}-{}", c[2], c[3], c[4])))
                .set("mb", r1(length / 1e6, 1));
            if inside != 0.0 {
                row.put("masked_mb", r1(inside / 1e6, 1));
            }
            big.push(row.val());
        }
    }
    let long_bp: f64 = per_chrom.values().fold(0.0, |a, v| a + v.0);
    let f_roh = long_bp / autosome_unmasked_bp() as f64;
    let mut upd = Value::Null;
    // The first chromosome of the most bases, as max() picks.
    let top =
        per_chrom.iter().fold(
            None::<(&String, &(f64, Vec<(i128, i128)>))>,
            |best, kv| match best {
                Some(b) if b.1.0 >= kv.1.0 => Some(b),
                _ => Some(kv),
            },
        );
    if let Some((chrom, (bp, runs))) = top
        && *bp >= UPD_MIN_BP
        && *bp > UPD_SHARE * long_bp
    {
        upd = Obj::new()
            .set("chrom", s(chrom))
            .set("mb", r1(bp / 1e6, 1))
            .set("share", r1(bp / long_bp, 2))
            .set("imprinted", imprinted(chrom, runs))
            .val();
    }
    let relationship = if upd == Value::Null {
        s(roh_band(f_roh))
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
        .set("f_roh", r1(f_roh, 4))
        .set("relationship", relationship)
        .set("upd_flag", upd);
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
