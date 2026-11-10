//! `og pgx cyp2d6-bed|cyp2d6-measure|cyp2d6-check`: is the read depth at
//! CYP2D6 fit for a copy-number call? The port of `bin/cyp2d6_depth_check.py`.
//!
//! pypgx (step 32) and Cyrius (step 21) call CYP2D6 copy number from read
//! depth. Before either calls, mosdepth measures the mean depth over CYP2D6
//! and two 50 kb flanks outside the CYP2D6-CYP2D8 cluster, once for all
//! reads and once for reads with MAPQ >= 1 (callers ignore MAPQ 0).
//!
//! - `cyp2d6-bed` prints those regions (GRCh38, 0-based BED, name in column 4).
//! - `cyp2d6-measure --bam BAM --out SHARES.tsv` counts the depth over
//!   Cyrius's CYP2D6 region (with REP6) and over the gene body alone, for all
//!   reads, MAPQ >= 1 and MAPQ >= 20, from the BAM itself.
//! - `cyp2d6-check --all Q0 --mapq1 Q1 [--measure SHARES.tsv] --out CHECK.tsv`
//!   reads the two mosdepth region files (plain or gzip) and writes the
//!   verdict, `ok` or `unreliable`; the shares never change it. Exit 0
//!   whatever the verdict, 2 on unreadable input.
//!
//! One deliberate change from the Python: a region listed twice in a
//! mosdepth file counts once, and listed twice with two different depths it
//! is refused (exit 2). The Python added every copy to the totals, so two
//! depths for one region were silently averaged.

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::io::{self, Write};
use std::path::Path;

use og_io::{AtomicFile, Lines, Utf8};
use pycompat::fmt::fixed;
use pycompat::num::{parse_float, parse_int};
use pycompat::text::strip;

use crate::args::Args;
use crate::{Bam, blocks};

const CHROM: &str = "chr22";

/// The kind of a region: CYP2D6 as Cyrius 1.1.1 defines it (with REP6), or
/// one of the flanks `scripts/ci/alt-depth-ab.sh` measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Kind {
    Gene,
    Flank,
}

/// GRCh38, 0-based, half-open, in the order `cyp2d6-bed` prints them.
const REGIONS: [(&str, i64, i64, &str); 3] = [
    (CHROM, 42_050_000, 42_100_000, "flank"),
    (CHROM, 42_123_192, 42_132_032, "CYP2D6"),
    (CHROM, 42_200_000, 42_250_000, "flank"),
];
const GENE: (i64, i64) = (42_123_192, 42_132_032);

/// The gene body without REP6: CYP2D6's MANE Select transcript
/// (ENST00000645361, NM_000106.6), GRCh38 chr22:42,126,499-42,130,810.
const BODY: (i64, i64) = (42_126_498, 42_130_810);

/// The MAPQ floors `measure` counts at.
const MEASURE_MAPQ: [u8; 3] = [0, 1, 20];

/// The reads mosdepth leaves out by default (its -F 1796): unmapped,
/// secondary, QC-failed and duplicate. `measure` leaves out the same.
const SKIP_FLAGS: u16 = 0x4 | 0x100 | 0x200 | 0x400;

const NOT_MEASURED: &str = "not measured";
const OK: &str = "CYP2D6 depth is fit for a copy-number call";
const UNRELIABLE: &str = "CYP2D6 copy number unreliable: reads are multi-mapped \
                          (was this BAM aligned to a reference with ALT contigs?)";
const NO_DEPTH: &str = "CYP2D6 copy number unreliable: no reads in the flanks of CYP2D6 \
                        (chr22:42.05-42.25 Mb), so its depth cannot be compared";

/// A depth fit for a copy-number call has at least this share of its reads
/// at MAPQ >= 1, at CYP2D6 and in the flanks.
const MIN_MAPQ1_SHARE: f64 = 0.6;
/// Below this ratio of CYP2D6 to flank depth (all reads) the gene has too
/// few reads for its MAPQ share to mean anything: both copies deleted.
const MIN_GENE_RATIO: f64 = 0.15;

/// The verdict on the four mean depths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    Ok,
    /// More than 40% of the reads at CYP2D6 or in its flanks have MAPQ 0.
    MultiMapped,
    /// No reads in the flanks, so nothing to compare CYP2D6 with.
    NoFlankDepth,
}

impl Verdict {
    pub fn status(self) -> &'static str {
        match self {
            Verdict::Ok => "ok",
            _ => "unreliable",
        }
    }

    pub fn message(self) -> &'static str {
        match self {
            Verdict::Ok => OK,
            Verdict::MultiMapped => UNRELIABLE,
            Verdict::NoFlankDepth => NO_DEPTH,
        }
    }
}

/// Mean depths over CYP2D6 and over its flanks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Depths {
    pub gene: f64,
    pub flank: f64,
}

/// The verdict for the depths of all reads and of MAPQ >= 1.
///
/// Unreliable when more than 40% of the reads at CYP2D6 or in its flanks have
/// MAPQ 0: the aligner placed them on more than one sequence. A deletion
/// leaves few reads at CYP2D6, but those map uniquely and the flanks are
/// untouched, so it passes; below a gene/flank ratio of 0.15 the gene share
/// is not read at all.
pub fn assess(all: Depths, q1: Depths) -> Verdict {
    if all.flank <= 0.0 {
        return Verdict::NoFlankDepth;
    }
    let ratio_all = all.gene / all.flank;
    let share_gene = if all.gene > 0.0 {
        q1.gene / all.gene
    } else {
        1.0
    };
    let share_flank = q1.flank / all.flank;
    if share_flank < MIN_MAPQ1_SHARE
        || (ratio_all >= MIN_GENE_RATIO && share_gene < MIN_MAPQ1_SHARE)
    {
        Verdict::MultiMapped
    } else {
        Verdict::Ok
    }
}

fn lines(path: &str) -> Result<Lines<Box<dyn io::BufRead>>, String> {
    let f = og_io::open(path).map_err(|e| format!("{path}: {e}"))?;
    Ok(Lines::new(f, Utf8::Strict))
}

fn int(s: &str, path: &str) -> Result<i64, String> {
    parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| format!("{path}: not an integer: {s:?}"))
}

/// The length-weighted mean depth of CYP2D6 and of the flanks from a
/// mosdepth regions file (chrom, start, end, name, mean). Rows named neither
/// are skipped; the regions must be the ones `cyp2d6-bed` prints. A row
/// listed twice counts once; the same region with two depths is an error.
pub fn read_regions(path: &str) -> Result<Depths, String> {
    // Each region with the bits of its mean depth.
    let mut seen: HashMap<(String, i64, i64, Kind), u64> = HashMap::new();
    let (mut total, mut length) = ([0.0f64; 2], [0i64; 2]);
    for line in lines(path)? {
        let line = line.map_err(|e| format!("{path}: {e}"))?;
        if strip(&line, None).is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            return Err(format!(
                "{path}: expected chrom, start, end, name, mean; got {:?}",
                strip(&line, None)
            ));
        }
        let kind = if f[3] == "CYP2D6" {
            Kind::Gene
        } else if pycompat::text::lower(f[3]).contains("flank") {
            Kind::Flank
        } else {
            continue;
        };
        let (start, end) = (int(f[1], path)?, int(f[2], path)?);
        let mean = parse_float(f[4]).map_err(|_| format!("{path}: not a number: {:?}", f[4]))?;
        match seen.entry((f[0].to_owned(), start, end, kind)) {
            Entry::Occupied(o) if *o.get() == mean.to_bits() => continue,
            Entry::Occupied(_) => {
                return Err(format!(
                    "{path}: {}:{start}-{end} ({}) is listed twice with two depths",
                    f[0], f[3]
                ));
            }
            Entry::Vacant(v) => {
                v.insert(mean.to_bits());
            }
        }
        let k = kind as usize;
        total[k] += mean * (end - start) as f64;
        length[k] += end - start;
    }
    let have = |k: Kind| seen.keys().any(|s| s.3 == k);
    let missing: Vec<&str> = [("flank", Kind::Flank), ("gene", Kind::Gene)]
        .into_iter()
        .filter(|&(_, k)| !have(k))
        .map(|(n, _)| n)
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "{path}: no {} region (name CYP2D6, or one containing 'flank')",
            missing.join(" or ")
        ));
    }
    let want: HashSet<(String, i64, i64, Kind)> = REGIONS
        .iter()
        .map(|&(c, s, e, n)| {
            let k = if n == "CYP2D6" {
                Kind::Gene
            } else {
                Kind::Flank
            };
            (c.to_owned(), s, e, k)
        })
        .collect();
    if seen.keys().cloned().collect::<HashSet<_>>() != want {
        return Err(format!(
            "{path}: its regions differ from the ones `cyp2d6_depth_check.py bed` prints"
        ));
    }
    let mean = |k: Kind| total[k as usize] / length[k as usize] as f64;
    Ok(Depths {
        gene: mean(Kind::Gene),
        flank: mean(Kind::Flank),
    })
}

/// One read as `measure` counts it: its flag, MAPQ and aligned blocks.
pub struct Counted {
    pub flag: u16,
    pub mapq: u8,
    pub blocks: Vec<(i64, i64)>,
}

/// Mean depth over `[start, end)` of the aligned bases of the reads at
/// `min_mapq` or more, without the reads mosdepth skips by default. Both
/// mates of an overlapping pair count, as in the Python: these depths only
/// make shares, which do not feed the verdict.
pub fn mean_depth(reads: &[Counted], start: i64, end: i64, min_mapq: u8) -> f64 {
    let total: i64 = reads
        .iter()
        .filter(|r| r.flag & SKIP_FLAGS == 0 && r.mapq >= min_mapq)
        .flat_map(|r| &r.blocks)
        .map(|&(s, e)| (e.min(end) - s.max(start)).max(0))
        .sum();
    total as f64 / (end - start) as f64
}

/// `(region, min_mapq, mean depth)` for Cyrius's region and the gene body.
pub fn measure_rows(reads: &[Counted]) -> Vec<(&'static str, u8, f64)> {
    [("gene", GENE), ("body", BODY)]
        .into_iter()
        .flat_map(|(name, (s, e))| {
            MEASURE_MAPQ
                .into_iter()
                .map(move |q| (name, q, mean_depth(reads, s, e, q)))
        })
        .collect()
}

fn measure(bam: &str, out: &str) -> Result<(), String> {
    let mut b = Bam::open(Path::new(bam)).map_err(|e| format!("{bam}: {e}"))?;
    let chrom = [CHROM, &CHROM[3..]]
        .into_iter()
        .find(|c| b.tid(c).is_some())
        .ok_or_else(|| format!("{bam} has no {CHROM}"))?;
    let mut reads = Vec::new();
    for r in b
        .fetch(chrom, GENE.0, GENE.1)
        .map_err(|e| format!("{bam}: {e}"))?
    {
        let r = r.map_err(|e| format!("{bam}: {e}"))?;
        reads.push(Counted {
            flag: r.flags(),
            mapq: r.mapq(),
            blocks: blocks(&r),
        });
    }
    let mut text = String::from("region\tmin_mapq\tmean_depth\n");
    for (n, q, d) in measure_rows(&reads) {
        text.push_str(&format!("{n}\t{q}\t{}\n", fixed(d, 4)));
    }
    write_atomic(out, &text)
}

fn write_atomic(path: &str, text: &str) -> Result<(), String> {
    let mut f = AtomicFile::create(path).map_err(|e| format!("cannot write {path}: {e}"))?;
    f.write_all(text.as_bytes())
        .and_then(|()| f.commit())
        .map_err(|e| format!("cannot write {path}: {e}"))
}

/// `((region, min_mapq), mean depth)` rows of `measure`'s file.
type Measured = Vec<((String, i64), f64)>;

/// The rows of `measure`'s file; a repeated row replaces the earlier one.
fn read_measure(path: &str) -> Result<Measured, String> {
    let mut it = lines(path)?;
    let head = it
        .next()
        .transpose()
        .map_err(|e| format!("{path}: {e}"))?
        .unwrap_or_default();
    if head != "region\tmin_mapq\tmean_depth" {
        return Err(format!(
            "{path}: not a file of `cyp2d6_depth_check.py measure`"
        ));
    }
    let mut out: Measured = Vec::new();
    for line in it {
        let line = line.map_err(|e| format!("{path}: {e}"))?;
        if strip(&line, None).is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let [n, q, d] = f[..] else {
            return Err(format!(
                "{path}: expected region, min_mapq, mean_depth: {line:?}"
            ));
        };
        let key = (n.to_owned(), int(q, path)?);
        let d = parse_float(d).map_err(|_| format!("{path}: not a number: {d:?}"))?;
        match out.iter_mut().find(|(k, _)| *k == key) {
            Some(slot) => slot.1 = d,
            None => out.push((key, d)),
        }
    }
    Ok(out)
}

/// `num/den` with two decimals, `no reads` when `den` is 0, `not measured`
/// when either is missing.
fn share(num: Option<f64>, den: Option<f64>) -> String {
    match (num, den) {
        (Some(n), Some(d)) if d > 0.0 => fixed(n / d, 2),
        (Some(_), Some(_)) => "no reads".to_owned(),
        _ => NOT_MEASURED.to_owned(),
    }
}

fn check(all: &str, mapq1: &str, measured: Option<&str>, out: &str) -> Result<String, String> {
    let d_all = read_regions(all)?;
    let d_q1 = read_regions(mapq1)?;
    let m = measured.map(read_measure).transpose()?.unwrap_or_default();
    let get = |n: &str, q: i64| {
        m.iter()
            .find(|((k, kq), _)| k == n && *kq == q)
            .map(|x| x.1)
    };
    let v = assess(d_all, d_q1);
    let body_all = get("body", 0);
    let rows: [(&str, String); 11] = [
        ("status", v.status().to_owned()),
        ("message", v.message().to_owned()),
        ("gene_depth_all", fixed(d_all.gene, 2)),
        ("gene_depth_mapq1", fixed(d_q1.gene, 2)),
        ("flank_depth_all", fixed(d_all.flank, 2)),
        ("flank_depth_mapq1", fixed(d_q1.flank, 2)),
        // The share the verdict reads, then what `measure` counted.
        ("gene_mapq1_share", share(Some(d_q1.gene), Some(d_all.gene))),
        ("gene_mapq20_share", share(get("gene", 20), get("gene", 0))),
        (
            "body_depth_all",
            body_all.map_or_else(|| NOT_MEASURED.to_owned(), |b| fixed(b, 2)),
        ),
        ("body_mapq1_share", share(get("body", 1), body_all)),
        ("body_mapq20_share", share(get("body", 20), body_all)),
    ];
    let mut text = String::from("metric\tvalue\n");
    for (k, val) in &rows {
        text.push_str(&format!("{k}\t{val}\n"));
    }
    write_atomic(out, &text)?;
    Ok(format!(
        "CYP2D6 depth check: {}: {}\n  mean depth, all reads: CYP2D6 {}, flanks {}; MAPQ >= 1: CYP2D6 {}, flanks {}\n  share of the depth at MAPQ >= 1 / >= 20: Cyrius's region (with REP6) {} / {}; gene body without REP6 {} / {}\n",
        v.status(),
        v.message(),
        rows[2].1,
        rows[4].1,
        rows[3].1,
        rows[5].1,
        rows[6].1,
        rows[7].1,
        rows[9].1,
        rows[10].1,
    ))
}

/// `og pgx cyp2d6-<verb> ARGS`: the exit code; usage errors are exit 2.
pub fn cli(verb: &str, argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let run = || -> Result<Result<String, String>, String> {
        Ok(match verb {
            "bed" => {
                Args::parse(argv, &[], &[])?;
                Ok(REGIONS
                    .iter()
                    .map(|(c, s, e, n)| format!("{c}\t{s}\t{e}\t{n}\n"))
                    .collect())
            }
            "measure" => {
                let a = Args::parse(argv, &["--bam", "--out"], &[])?;
                let (bam, o) = (a.required("--bam")?, a.required("--out")?);
                measure(bam, o).map(|()| String::new())
            }
            "check" => {
                let a = Args::parse(argv, &["--all", "--mapq1", "--measure", "--out"], &[])?;
                let (all, q1, o) = (
                    a.required("--all")?,
                    a.required("--mapq1")?,
                    a.required("--out")?,
                );
                check(all, q1, a.get("--measure"), o)
            }
            other => return Err(format!("unknown verb 'cyp2d6-{other}'")),
        })
    };
    match run() {
        Err(usage) => {
            let _ = writeln!(
                err,
                "usage: og pgx cyp2d6-bed | cyp2d6-measure --bam BAM --out OUT | \
                 cyp2d6-check --all Q0 --mapq1 Q1 [--measure SHARES] --out OUT\nog: error: {usage}"
            );
            2
        }
        Ok(Err(e)) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
        Ok(Ok(text)) => match out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "ERROR: {e}");
                2
            }
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(gene: f64, flank: f64) -> Depths {
        Depths { gene, flank }
    }

    #[test]
    fn the_verdicts() {
        assert_eq!(assess(d(10.0, 0.0), d(10.0, 0.0)), Verdict::NoFlankDepth);
        assert_eq!(assess(d(30.0, 30.0), d(29.0, 29.0)), Verdict::Ok);
        // A third of the flank reads at MAPQ 0: multi-mapped.
        assert_eq!(assess(d(30.0, 30.0), d(29.0, 15.0)), Verdict::MultiMapped);
        // The gene share is read only above the gene ratio floor.
        assert_eq!(assess(d(3.0, 30.0), d(0.0, 29.0)), Verdict::Ok);
        assert_eq!(assess(d(6.0, 30.0), d(0.0, 29.0)), Verdict::MultiMapped);
        // No gene depth at all: the gene share is taken as 1.
        assert_eq!(assess(d(0.0, 30.0), d(0.0, 29.0)), Verdict::Ok);
    }

    #[test]
    fn shares() {
        assert_eq!(share(Some(1.0), Some(3.0)), "0.33");
        assert_eq!(share(Some(1.0), Some(0.0)), "no reads");
        assert_eq!(share(None, Some(1.0)), NOT_MEASURED);
    }

    #[test]
    fn depth_counts_aligned_bases_above_the_floor_only() {
        let reads = [
            Counted {
                flag: 0,
                mapq: 60,
                blocks: vec![(0, 10)],
            },
            Counted {
                flag: 0,
                mapq: 0,
                blocks: vec![(5, 15)],
            },
            Counted {
                flag: 0x400,
                mapq: 60,
                blocks: vec![(0, 10)],
            },
        ];
        assert_eq!(mean_depth(&reads, 0, 10, 0), 1.5);
        assert_eq!(mean_depth(&reads, 0, 10, 1), 1.0);
    }
}
