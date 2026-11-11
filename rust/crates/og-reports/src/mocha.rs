//! `og mocha report` (step 51): the mosaic chromosomal alterations MoChA
//! calls, with MoChA's own UK Biobank filters and this pipeline's floor.
//!
//! A mosaic chromosomal alteration (mCA) is a large gain, loss or
//! copy-neutral loss of heterozygosity carried by a share of the blood
//! cells. They rise with age and go with a higher risk of a later blood
//! cancer (Loh et al. 2018, PMID 29995854). MoChA calls them from a phased
//! VCF: in the cells that carry the event one parental haplotype is lost or
//! gained, so the allele balance leans the same way along the haplotype.
//!
//! ```text
//! og mocha report --sample S --calls CALLS.tsv --stats STATS.tsv --sites SITES.tsv \
//!     --phasing beagle|input --out-calls S_mocha.tsv --out-summary S_mocha_summary.tsv
//! ```
//!
//! A call is reported when it is on chr1 to chr22, MoChA's filter of the
//! UK Biobank callset keeps it (not a germline CNP; not likely
//! constitutional; not weak; not a likely germline duplication), and its
//! cell fraction is at or above the floor: the larger of 0.05 and
//! 6 / sqrt(D x H), D MoChA's median depth and H the heterozygous sites per
//! 10 Mb. The sample is not assessed when MoChA's sample filter fails
//! (call rate under 0.97, phased BAF autocorrelation over 0.03) or cannot
//! be applied. Exit 0 whatever the result; 2 on unreadable input.
//!
//! NaN is kept as Python has it: every comparison with NaN is false, so a
//! call with no BAF deviation, no N50 and no concordance LOD is dropped as
//! likely constitutional, as MoChA's own filter expression does.
//!
//! What changes on purpose: `calls_mocha` counts the calls of this sample
//! (the Python counted every row of a multi-sample table, other samples'
//! included; a pipeline run gives MoChA one sample, so its value does not
//! change), and a calls row that ends before a needed column, or a sites
//! value that is not ASCII digits, is an input error (exit 2), not a
//! traceback.

use crate::{Error, Outcome};
use og_filters::args::Spec;
use og_io::Utf8;
use og_tsv::{Options, Quoting};
use std::io::Write;

const FLOOR_MIN: f64 = 0.05;
const FLOOR_Z: f64 = 3.0;
const FLOOR_MB: f64 = 10.0;
const CALL_RATE_MIN: f64 = 0.97;
const BAF_AUTO_MAX: f64 = 0.03;

pub const CALL_COLUMNS: [&str; 14] = [
    "chrom",
    "start",
    "end",
    "length_mb",
    "arms",
    "type",
    "cell_fraction",
    "bdev",
    "rel_cov",
    "n_hets",
    "lod_baf_phase",
    "lod_lrr_baf",
    "reported",
    "reason",
];
const NEED_CALLS: [&str; 14] = [
    "chrom",
    "length",
    "p_arm",
    "q_arm",
    "n50_hets",
    "bdev",
    "bdev_se",
    "rel_cov",
    "lod_baf_phase",
    "lod_lrr_baf",
    "lod_baf_conc",
    "n_hets",
    "type",
    "cf",
];
const NEED_STATS: [&str; 5] = ["sample_id", "call_rate", "cov_median", "baf_auto", "n_hets"];

/// chr1 to chr22, with or without the prefix.
fn autosome(c: &str) -> bool {
    let n = c.strip_prefix("chr").unwrap_or(c);
    crate::digits(n) && !n.starts_with('0') && n.parse::<u8>().is_ok_and(|v| (1..=22).contains(&v))
}

/// MoChA's number: NaN for `nan`, `.`, an empty or a missing cell.
fn num(x: Option<&str>) -> f64 {
    x.and_then(crate::float).unwrap_or(f64::NAN)
}

/// `f"{x:.{d}f}"`, "." for NaN.
fn fmt(x: f64, d: usize) -> String {
    if x.is_nan() {
        ".".into()
    } else {
        crate::fixed(x, d)
    }
}

/// MoChA's call type, as the table writes it.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Kind {
    Loss,
    Gain,
    CnLoh,
    Undetermined,
    Other(String),
}

impl Kind {
    fn read(s: &str) -> Kind {
        match s {
            "Loss" => Kind::Loss,
            "Gain" => Kind::Gain,
            "CN-LOH" => Kind::CnLoh,
            "Undetermined" => Kind::Undetermined,
            o => Kind::Other(o.into()),
        }
    }

    fn word(&self) -> &str {
        match self {
            Kind::Loss => "loss",
            Kind::Gain => "gain",
            Kind::CnLoh => "CN-LOH",
            Kind::Undetermined => "undetermined",
            Kind::Other(s) => s,
        }
    }
}

/// One chromosome of the sites table.
#[derive(Debug, Clone)]
pub struct Sites {
    chrom: String,
    first: i128,
    last: i128,
    hets: i128,
    phased: i128,
}

/// A table of MoChA's, checked for the columns it needs.
fn read_table(path: &str, need: &[&str]) -> Result<Vec<og_tsv::Row>, Error> {
    let opts = Options::default()
        .quoting(Quoting::None)
        .utf8(Utf8::Replace);
    let (header, rows) = crate::read_table(path, opts)?;
    let missing: Vec<&str> = need
        .iter()
        .copied()
        .filter(|k| !header.iter().any(|h| h == k))
        .collect();
    if !missing.is_empty() {
        return Err(Error(format!(
            "{path} lacks the column{} {}",
            if missing.len() > 1 { "s" } else { "" },
            missing.join(", ")
        )));
    }
    Ok(rows)
}

fn read_sites(path: &str) -> Result<Vec<Sites>, Error> {
    let mut out = Vec::new();
    for (n, line) in crate::read_lines(path, Utf8::Replace)?.iter().enumerate() {
        if crate::strip(line).is_empty() {
            continue;
        }
        let c: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(line)
            .split('\t')
            .collect();
        let v: Option<Vec<i128>> = (c.len() == 6)
            .then(|| {
                c[1..]
                    .iter()
                    .map(|x| crate::digits(x).then(|| crate::int(x)).flatten())
                    .collect()
            })
            .flatten();
        let Some(v) = v else {
            return Err(Error(format!(
                "{path} line {}: want chrom, first, last, sites, hets, phased hets",
                n + 1
            )));
        };
        out.push(Sites {
            chrom: c[0].into(),
            first: v[0],
            last: v[1],
            hets: v[3],
            phased: v[4],
        });
    }
    Ok(out)
}

/// The read-sampling floor: 6 / sqrt(D x H); NaN without data.
fn floor_read(depth: f64, het_per_10mb: f64) -> f64 {
    if !(depth > 0.0 && het_per_10mb > 0.0) {
        return f64::NAN;
    }
    4.0 * FLOOR_Z * 0.5 / (depth * het_per_10mb).sqrt()
}

/// `20p`, `20q` or `20p and 20q`, with the telomere it reaches.
fn arms(chrom: &str, p_arm: Option<&str>, q_arm: Option<&str>) -> String {
    let c = chrom.strip_prefix("chr").unwrap_or(chrom);
    let both = [("p", p_arm), ("q", q_arm)];
    let on: Vec<&str> = both
        .iter()
        .filter(|(_, v)| matches!(v, Some("Y" | "C" | "T")))
        .map(|(a, _)| *a)
        .collect();
    let tel: Vec<&str> = both
        .iter()
        .filter(|(_, v)| *v == Some("T"))
        .map(|(a, _)| *a)
        .collect();
    if on.is_empty() {
        return format!("{c} (interstitial)");
    }
    let mut text = on
        .iter()
        .map(|a| format!("{c}{a}"))
        .collect::<Vec<_>>()
        .join(" and ");
    if !tel.is_empty() {
        text.push_str(", to the ");
        text.push_str(
            &tel.iter()
                .map(|a| format!("{a} telomere"))
                .collect::<Vec<_>>()
                .join(" and "),
        );
    }
    text
}

/// The reason MoChA's UK Biobank filter drops a call, or `None`.
// The negated comparisons are on purpose: a NaN fails every comparison, as
// in the Python, so `!(x > 10.0)` holds for a NaN and `x <= 10.0` would not.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn mocha_filter(r: &og_tsv::Row) -> Option<String> {
    let t = r.get("type").unwrap_or("");
    if t.starts_with("CNP") {
        return Some(format!(
            "germline copy-number polymorphism (MoChA type {t})"
        ));
    }
    let n = |k| num(r.get(k));
    let (bdev, n50, conc) = (n("bdev"), n("n50_hets"), n("lod_baf_conc"));
    if !(bdev < 0.1 || n50 < 2e5 || (!conc.is_nan() && conc > 10.0)) {
        return Some(
            "likely constitutional (BAF deviation 0.1 or more with sparse heterozygous sites)"
                .into(),
        );
    }
    if n("bdev_se").is_nan() && !(n("lod_baf_phase") > 10.0) {
        return Some(
            "weak: no BAF deviation estimate and a BAF-and-phase LOD of 10 or less".into(),
        );
    }
    let (rel, length) = (n("rel_cov"), n("length"));
    if !(rel < 2.1
        || bdev < 0.05
        || (length > 5e5 && bdev < 0.1 && rel < 2.5)
        || (length > 5e6 && bdev < 0.15))
    {
        return Some("likely germline duplication (relative coverage 2.1 or more)".into());
    }
    None
}

/// The report: call rows (by [`CALL_COLUMNS`]) and the summary rows.
pub type Report = (Vec<Vec<String>>, Vec<(&'static str, String)>);

/// The report of `sample` from MoChA's rows, the sites and the phasing.
pub fn report(
    sample: &str,
    calls: &[og_tsv::Row],
    stats: &[og_tsv::Row],
    sites: &[Sites],
    phasing: &str,
) -> Result<Report, Error> {
    let st = stats
        .iter()
        .find(|r| r.get("sample_id") == Some(sample))
        .or(if stats.len() == 1 {
            stats.first()
        } else {
            None
        })
        .ok_or_else(|| Error(format!("MoChA's stats table has no row for {sample}")))?;
    let auto: Vec<&Sites> = sites.iter().filter(|s| autosome(&s.chrom)).collect();
    let span: i128 = auto.iter().map(|s| (s.last - s.first + 1).max(0)).sum();
    let hets: i128 = auto.iter().map(|s| s.hets).sum();
    let phased: i128 = auto.iter().map(|s| s.phased).sum();
    let depth = num(st.get("cov_median"));
    let h10 = if span != 0 {
        hets as f64 * FLOOR_MB * 1e6 / span as f64
    } else {
        f64::NAN
    };
    let fr = floor_read(depth, h10);
    let floor = if fr.is_nan() {
        f64::NAN
    } else if fr > FLOOR_MIN {
        fr
    } else {
        FLOOR_MIN
    };
    let (call_rate, baf_auto) = (num(st.get("call_rate")), num(st.get("baf_auto")));
    let why = if auto.is_empty() || hets == 0 {
        Some("no heterozygous site on chr1 to chr22 reached MoChA".to_owned())
    } else if phased == 0 {
        Some("no heterozygous site is phased, so MoChA cannot follow a haplotype".into())
    } else if floor.is_nan() {
        Some("MoChA gave no median depth, so no floor can be set".into())
    } else if call_rate.is_nan() || baf_auto.is_nan() {
        Some("MoChA gave no call rate or phased BAF autocorrelation, so its sample filter cannot be applied".into())
    } else if call_rate < CALL_RATE_MIN {
        Some(format!(
            "call rate {} is under {} (MoChA's sample filter)",
            crate::fixed(call_rate, 4),
            pycompat::fmt::repr(CALL_RATE_MIN)
        ))
    } else if baf_auto > BAF_AUTO_MAX {
        Some(format!(
            "phased BAF autocorrelation {} is over {} (MoChA's sample filter): noisy data",
            crate::fixed(baf_auto, 4),
            pycompat::fmt::repr(BAF_AUTO_MAX)
        ))
    } else {
        None
    };
    let mut rows = Vec::new();
    for r in calls {
        if stats.len() > 1 && r.get("sample_id").is_some_and(|s| s != sample) {
            continue;
        }
        for k in NEED_CALLS {
            r.text(k)?;
        }
        let chrom = r.get("chrom").unwrap_or("");
        let cf = num(r.get("cf"));
        let reason = if !autosome(chrom) {
            Some("not an autosome (mosaic loss of Y and X: step 33)".to_owned())
        } else if let Some(f) = mocha_filter(r) {
            Some(f)
        } else if cf.is_nan() {
            Some("cell fraction not estimated".into())
        } else if !floor.is_nan() && cf < floor {
            Some(format!(
                "cell fraction {} under the floor {}",
                crate::fixed(cf, 3),
                crate::fixed(floor, 3)
            ))
        } else if why.is_some() {
            Some("sample not assessed".into())
        } else {
            None
        };
        let prefixed = |p: &str| {
            r.dict()
                .iter()
                .find(|(k, _)| k.starts_with(p))
                .map_or(".".to_owned(), |(_, v)| v.unwrap_or("None").to_owned())
        };
        let raw = |k: &str| r.get(k).unwrap_or("None").to_owned();
        rows.push(vec![
            chrom.to_owned(),
            prefixed("beg_"),
            prefixed("end_"),
            fmt(num(r.get("length")) / 1e6, 2),
            arms(chrom, r.get("p_arm"), r.get("q_arm")),
            Kind::read(r.get("type").unwrap_or("")).word().to_owned(),
            fmt(cf, 3),
            fmt(num(r.get("bdev")), 4),
            fmt(num(r.get("rel_cov")), 4),
            raw("n_hets"),
            raw("lod_baf_phase"),
            raw("lod_lrr_baf"),
            if reason.is_some() { "no" } else { "yes" }.into(),
            reason.unwrap_or_else(|| ".".into()),
        ]);
    }
    let reported = rows.iter().filter(|r| r[12] == "yes").count();
    let summary = vec![
        (
            "status",
            if why.is_some() {
                "not assessed"
            } else {
                "assessed"
            }
            .to_owned(),
        ),
        ("message", why.unwrap_or_else(|| ".".into())),
        ("depth", fmt(depth, 1)),
        ("het_sites", hets.to_string()),
        ("phased_het_sites", phased.to_string()),
        ("het_per_10mb", fmt(h10, 0)),
        ("analysed_mb", fmt(span as f64 / 1e6, 1)),
        (
            "chromosomes",
            if auto.is_empty() {
                ".".into()
            } else {
                auto.iter()
                    .map(|s| s.chrom.as_str())
                    .collect::<Vec<_>>()
                    .join(",")
            },
        ),
        ("call_rate", fmt(call_rate, 4)),
        ("baf_auto", fmt(baf_auto, 4)),
        ("floor_read", fmt(fr, 3)),
        ("floor_min", fmt(FLOOR_MIN, 2)),
        ("floor", fmt(floor, 3)),
        ("calls_mocha", rows.len().to_string()),
        ("calls_reported", reported.to_string()),
        ("phasing", phasing.to_owned()),
    ];
    Ok((rows, summary))
}

const SPEC: Spec = Spec {
    prog: "og mocha report",
    usage: "usage: og mocha report [-h] --sample SAMPLE --calls CALLS --stats STATS --sites SITES --phasing {beagle,input} --out-calls OUT_CALLS --out-summary OUT_SUMMARY\n",
    options: &[
        "--sample",
        "--calls",
        "--stats",
        "--sites",
        "--phasing",
        "--out-calls",
        "--out-summary",
    ],
    positionals: &[],
    required: &[
        "--sample",
        "--calls",
        "--stats",
        "--sites",
        "--phasing",
        "--out-calls",
        "--out-summary",
    ],
};

/// `og mocha VERB ARGS`: the exit code.
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if verb != Some("report") {
        return crate::unknown_verb("mocha", "report", verb, err);
    }
    crate::run(&SPEC, rest, out, err, |a| {
        let phasing = crate::opt(a, "--phasing");
        if phasing != "beagle" && phasing != "input" {
            return Err(og_filters::args::ArgError::Usage(format!(
                "argument --phasing: invalid choice: {} (choose from 'beagle', 'input')",
                og_filters::py_repr(phasing)
            )));
        }
        Ok(run(
            crate::opt(a, "--sample"),
            crate::opt(a, "--calls"),
            crate::opt(a, "--stats"),
            crate::opt(a, "--sites"),
            phasing,
            crate::opt(a, "--out-calls"),
            crate::opt(a, "--out-summary"),
        ))
    })
}

fn run(
    sample: &str,
    calls: &str,
    stats: &str,
    sites: &str,
    phasing: &str,
    out_calls: &str,
    out_summary: &str,
) -> Outcome {
    let bad = |e: Error| (e, 2);
    let c = read_table(calls, &NEED_CALLS).map_err(bad)?;
    let s = read_table(stats, &NEED_STATS).map_err(bad)?;
    let si = read_sites(sites).map_err(bad)?;
    let (rows, summary) = report(sample, &c, &s, &si, phasing).map_err(bad)?;
    let mut body = CALL_COLUMNS.join("\t") + "\n";
    for r in &rows {
        body.push_str(&r.join("\t"));
        body.push('\n');
    }
    crate::write_atomic(out_calls, body.as_bytes()).map_err(|e| (e, 1))?;
    let mut sum = String::from("metric\tvalue\n");
    for (k, v) in &summary {
        sum.push_str(&format!("{k}\t{v}\n"));
    }
    crate::write_atomic(out_summary, sum.as_bytes()).map_err(|e| (e, 1))?;
    let get = |k: &str| {
        summary
            .iter()
            .find(|(n, _)| *n == k)
            .map_or("", |(_, v)| v.as_str())
    };
    let msg = get("message");
    Ok(format!(
        "{sample}: {}; {} of {} MoChA calls reported (floor {}, depth {}x){}\n",
        get("status"),
        get("calls_reported"),
        get("calls_mocha"),
        get("floor"),
        get("depth"),
        if msg == "." {
            String::new()
        } else {
            format!(": {msg}")
        }
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arms_and_autosomes() {
        assert_eq!(
            arms("chr20", Some("N"), Some("T")),
            "20q, to the q telomere"
        );
        assert_eq!(arms("7", Some("C"), Some("Y")), "7p and 7q");
        assert_eq!(arms("chr1", Some("N"), Some("N")), "1 (interstitial)");
        assert!(autosome("chr22") && autosome("1") && !autosome("chr23") && !autosome("chrX"));
        assert!(!autosome("chr01") && !autosome("chr0"));
    }

    #[test]
    fn floor() {
        // 30x, 120k hets per 100 Mb: 6 / sqrt(30 * 12000) = 0.01, under 0.05.
        assert!((floor_read(30.0, 12_000.0) - 0.01).abs() < 1e-12);
        assert!(floor_read(f64::NAN, 1.0).is_nan());
        assert!(floor_read(0.0, 1.0).is_nan());
    }

    fn table(d: &crate::testutil::Dir, name: &str, text: &str, need: &[&str]) -> Vec<og_tsv::Row> {
        read_table(&d.put(name, text), need).unwrap()
    }

    const CALL_HEAD: &str = "sample_id\tchrom\tbeg_GRCh38\tend_GRCh38\tlength\tp_arm\tq_arm\tn_hets\tn50_hets\tbdev\tbdev_se\trel_cov\tlod_baf_phase\tlod_lrr_baf\tlod_baf_conc\ttype\tcf\n";
    const LOSS: &str = "S\tchr20\t30100000\t42300000\t12200000\tN\tY\t4210\t350000\t0.031\t0.004\t1.86\t48.2\t35.1\tnan\tLoss\t0.120\n";

    fn run_with(call_rate: &str, calls: &str) -> Report {
        let d = crate::testutil::Dir::new("mocha");
        let stats = format!(
            "sample_id\tcall_rate\tcov_median\tbaf_auto\tn_hets\nS\t{call_rate}\t31\t0.0124\t2100000\n"
        );
        let c = table(&d, "c.tsv", &format!("{CALL_HEAD}{calls}"), &NEED_CALLS);
        let s = table(&d, "s.tsv", &stats, &NEED_STATS);
        let sites: Vec<Sites> = (1..=22)
            .map(|i| Sites {
                chrom: format!("chr{i}"),
                first: 10001,
                last: 100_010_001,
                hets: 120_000,
                phased: 94_000,
            })
            .collect();
        report("S", &c, &s, &sites, "beagle").unwrap()
    }

    #[test]
    fn a_loss_is_reported() {
        let (rows, summary) = run_with("0.9981", LOSS);
        assert_eq!(rows[0][4], "20q");
        assert_eq!(rows[0][12], "yes", "{rows:?}");
        assert_eq!(summary[0].1, "assessed");
        assert_eq!(summary[12], ("floor", "0.050".to_owned()));
    }

    /// Red first: a call rate of exactly 0.97 meets MoChA's sample filter.
    #[test]
    fn red_first_a_call_rate_of_0_97_is_assessed() {
        let (rows, summary) = run_with("0.97", LOSS);
        assert_eq!(summary[0].1, "assessed");
        assert_eq!(rows[0][12], "yes");
        let (rows, summary) = run_with("0.9699", LOSS);
        assert_eq!(summary[0].1, "not assessed");
        assert_eq!(rows[0][13], "sample not assessed");
    }

    #[test]
    fn a_short_row_is_an_error() {
        let d = crate::testutil::Dir::new("mocha-short");
        let stats = table(
            &d,
            "s.tsv",
            "sample_id\tcall_rate\tcov_median\tbaf_auto\tn_hets\nS\t1\t31\t0\t1\n",
            &NEED_STATS,
        );
        let c = table(
            &d,
            "c.tsv",
            &format!("{CALL_HEAD}S\tchr20\t1\n"),
            &NEED_CALLS,
        );
        let e = report("S", &c, &stats, &[], "beagle").unwrap_err();
        assert!(
            e.0.ends_with("c.tsv:2: the row ends before column length"),
            "{e}"
        );
    }
}
