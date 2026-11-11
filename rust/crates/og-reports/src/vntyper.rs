//! `og vntyper report` (step 52): VNtyper's MUC1 VNTR result in plain words.
//!
//! ADTKD-MUC1 is caused by frameshifts inside the MUC1 coding VNTR, a repeat
//! of 60 bp units, most often one extra C in a run of seven (dupC; Kirby et
//! al., Nat Genet 2013, PMID 23396133). Short reads cannot be placed inside
//! the repeat, so standard callers do not see it. VNtyper (Saei et al.,
//! iScience 2023, PMID 37456840) genotypes it with k-mers (Kestrel); this
//! reads the run folder VNtyper writes and says positive, possible, negative
//! or not assessed, with the share of known positives VNtyper found at the
//! depth this sample has over the VNTR.
//!
//! ```text
//! og vntyper report --run-dir DIR --out TABLE [--exit-code N --log FILE]
//! ```
//!
//! DIR is VNtyper's `--output-dir` (`kestrel/kestrel_result.tsv` and
//! `coverage/coverage_summary.tsv`). `--exit-code` is VNtyper's own exit
//! status; a non-zero one, a missing file, or a Kestrel file with no result
//! row gives "not assessed". `--log`, VNtyper's console output, tells a BAM
//! with no read over MUC1 (VNtyper 2.0.47 stops on it with "FASTQ layout
//! 'empty'") from other errors. TABLE gets `metric<TAB>value` rows.
//!
//! What changes on purpose: a negative mean depth in VNtyper's coverage
//! table is an input error (exit 2); the Python stopped on an assertion.

use crate::{Error, Outcome};
use og_filters::args::Spec;
use og_tsv::{KvTable, Options, Utf8};
use std::io::Write;
use std::path::Path;

/// The rows of `FIELDS`, in order.
pub const FIELDS: [&str; 13] = [
    "status",
    "confidence",
    "variant",
    "flag",
    "alt_depth",
    "region_depth",
    "depth_score",
    "vntr_mean_depth",
    "vntr_uncovered_pct",
    "coverage_qc",
    "sensitivity",
    "message",
    "advice",
];

/// What VNtyper 2.0.47 says about sensitivity at depth (docs/user-guide/
/// configuration.md, "Coverage thresholds"): confirmed positive GRCh38
/// exomes downsampled and re-run, 42 samples at nine fractions each. Rows:
/// lowest mean depth over the VNTR, detected runs, runs.
const SENSITIVITY: [(f64, u32, u32); 3] = [(100.0, 160, 161), (50.0, 95, 108), (0.0, 48, 109)];
const SENSITIVITY_SOURCE: &str = "VNtyper 2.0.47's downsampling of 42 confirmed positive GRCh38 exomes (repeated runs, not independent tests)";

const ADVICE_CALL: &str = "A research call, not a diagnosis. VNtyper's authors ask for confirmation by another method: SNaPshot for dupC, long-read sequencing for other variants. ADTKD-MUC1 is dominant: discuss it with a nephrologist or a clinical geneticist.";
const ADVICE_NONE: &str = "A negative here lowers the chance of the frameshifts VNtyper looks for; it does not rule ADTKD-MUC1 out. A family history of kidney failure with gout or a bland urine sediment still deserves a clinical test (SNaPshot for dupC).";

const NO_READS: &str = "FASTQ layout 'empty'";
const REGION: &str = "chr1:155,184,000-155,194,000";
const KESTREL: &str = "kestrel/kestrel_result.tsv";
const COVERAGE: &str = "coverage/coverage_summary.tsv";

/// Kestrel's confidence of its best call.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Confidence {
    HighPrecisionStar,
    HighPrecision,
    LowPrecision,
    /// Anything else, `Negative` when Kestrel left it empty: no call.
    Other(String),
}

impl Confidence {
    fn read(s: &str) -> Confidence {
        match s {
            "High_Precision*" => Confidence::HighPrecisionStar,
            "High_Precision" => Confidence::HighPrecision,
            "Low_Precision" => Confidence::LowPrecision,
            "" => Confidence::Other("Negative".into()),
            other => Confidence::Other(other.into()),
        }
    }

    fn text(&self) -> &str {
        match self {
            Confidence::HighPrecisionStar => "High_Precision*",
            Confidence::HighPrecision => "High_Precision",
            Confidence::LowPrecision => "Low_Precision",
            Confidence::Other(s) => s,
        }
    }
}

/// Why a caveat makes a call "possible", not "positive".
#[derive(Debug, Clone, PartialEq)]
enum Caveat {
    Flagged(String),
    LowPrecision,
    Uncovered(f64),
    NoUncoveredMeasure,
}

impl Caveat {
    fn text(&self) -> String {
        match self {
            Caveat::Flagged(f) => format!("flagged: {f}"),
            Caveat::LowPrecision => "low precision".into(),
            Caveat::Uncovered(u) => format!("{}% of the VNTR has no read", crate::fixed(*u, 0)),
            Caveat::NoUncoveredMeasure => {
                "VNtyper measured no share of the VNTR without reads".into()
            }
        }
    }
}

/// The table, one value per field of [`FIELDS`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reading(Vec<(&'static str, String)>);

impl Reading {
    fn new() -> Reading {
        Reading(FIELDS.iter().map(|f| (*f, ".".to_owned())).collect())
    }

    fn set(&mut self, key: &str, value: impl Into<String>) {
        if let Some(slot) = self.0.iter_mut().find(|(k, _)| *k == key) {
            slot.1 = value.into();
        }
    }

    /// The value of `key`.
    pub fn get(&self, key: &str) -> &str {
        self.0
            .iter()
            .find(|(k, _)| *k == key)
            .map_or(".", |(_, v)| v)
    }

    fn not_assessed(mut self, message: String) -> Reading {
        self.set("status", "not assessed");
        self.set("message", message);
        self.set("advice", ADVICE_NONE);
        self
    }
}

/// A number of VNtyper's tables: `None` for a missing cell, text that is no
/// number, or NaN.
fn number(v: Option<&str>) -> Option<f64> {
    v.and_then(crate::float).filter(|x| !x.is_nan())
}

/// The sentence on VNtyper's sensitivity at the sample's mean depth.
fn sensitivity_line(mean: Option<f64>) -> Result<String, Error> {
    let Some(mean) = mean else {
        return Ok("Sensitivity not known: VNtyper measured no depth over the VNTR.".into());
    };
    let Some(&(low, hit, runs)) = SENSITIVITY.iter().find(|(low, ..)| mean >= *low) else {
        return Err(Error(format!(
            "VNtyper's coverage table gives a negative mean depth over the VNTR ({})",
            pycompat::fmt::repr(mean)
        )));
    };
    let band = if low == 100.0 {
        format!("at {}x or more", crate::fixed(low, 0))
    } else if low != 0.0 {
        format!("between {}x and 100x", crate::fixed(low, 0))
    } else {
        "below 50x".to_owned()
    };
    Ok(format!(
        "Mean depth over the MUC1 VNTR: {}x. At that depth ({band}) VNtyper found the known variant in {hit} of {runs} runs ({}%), in {SENSITIVITY_SOURCE}.",
        crate::fixed(mean, 0),
        crate::fixed(100.0 * f64::from(hit) / f64::from(runs), 0)
    ))
}

/// The data rows of a VNtyper TSV: its `#` banner lines and blank lines
/// skipped, the first line left the header.
fn rows(path: &str) -> Result<Vec<og_tsv::Row>, Error> {
    let opts = Options::default()
        .utf8(Utf8::Strict)
        .skip(|l| l.starts_with('#') || crate::strip(l).is_empty());
    Ok(crate::read_table(path, opts)?.1)
}

/// The reading of VNtyper's run folder `run_dir`.
pub fn read(run_dir: &str, exit_code: i64, log_text: &str) -> Result<Reading, Error> {
    let out = Reading::new();
    if exit_code != 0 {
        let why = if log_text.contains(NO_READS) {
            format!(
                "The BAM has no read in the MUC1 region ({REGION}), so VNtyper had nothing to read."
            )
        } else {
            format!("VNtyper stopped with exit code {exit_code}: see its log.")
        };
        return Ok(out.not_assessed(why));
    }
    let dir = Path::new(run_dir);
    let missing: Vec<&str> = [KESTREL, COVERAGE]
        .into_iter()
        .filter(|f| !dir.join(f).is_file())
        .collect();
    if !missing.is_empty() {
        return Ok(out.not_assessed(format!("VNtyper wrote no {}.", missing.join(" and no "))));
    }
    let path = |f: &str| dir.join(f).display().to_string();
    let cov = rows(&path(COVERAGE))?;
    let cov = cov.first();
    let cell = |r: Option<&og_tsv::Row>, k: &str| r.and_then(|r| r.get(k)).map(str::to_owned);
    let mean = number(cell(cov, "mean").as_deref());
    let uncovered = number(cell(cov, "percent_uncovered").as_deref());
    let qc = cell(cov, "coverage_qc")
        .map(|s| crate::strip(&s).to_owned())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "NOT_EVALUATED".into());
    let k = rows(&path(KESTREL))?;
    let Some(best) = k.first() else {
        return Ok(
            out.not_assessed("VNtyper wrote kestrel/kestrel_result.tsv with no result row.".into())
        );
    };
    let get = |c: &str| best.get(c).unwrap_or("");
    let conf = Confidence::read(crate::strip(get("Confidence")));
    let flag = crate::strip(get("Flag")).to_owned();
    let flagged = !flag.is_empty() && flag != "Not flagged";
    let covered = qc == "PASS" || qc == "REDUCED";
    let mut out = out;
    out.set("confidence", conf.text());
    out.set("coverage_qc", qc.clone());
    out.set("sensitivity", sensitivity_line(mean)?);
    out.set(
        "vntr_mean_depth",
        mean.map_or(".".into(), |m| crate::fixed(m, 2)),
    );
    out.set(
        "vntr_uncovered_pct",
        uncovered.map_or(".".into(), |u| crate::fixed(u, 2)),
    );
    if !matches!(conf, Confidence::Other(_)) {
        let name = crate::strip(get("Nomenclature"));
        let variant = if name.is_empty() {
            let v = crate::strip(best.get("Variant").filter(|v| !v.is_empty()).unwrap_or("."));
            if v.is_empty() { "." } else { v }
        } else {
            name
        };
        let or_dot = |c: &str| {
            let v = get(c);
            if v.is_empty() {
                ".".to_owned()
            } else {
                v.to_owned()
            }
        };
        out.set("variant", variant);
        out.set("flag", if flag.is_empty() { "." } else { &flag });
        out.set("alt_depth", or_dot("Estimated_Depth_AlternateVariant"));
        out.set(
            "region_depth",
            or_dot("Estimated_Depth_Variant_ActiveRegion"),
        );
        out.set("depth_score", or_dot("Depth_Score"));
        out.set("advice", ADVICE_CALL);
        let what = format!(
            "a MUC1 VNTR frameshift ({variant}), Kestrel confidence {}",
            conf.text()
        );
        let high = matches!(
            conf,
            Confidence::HighPrecision | Confidence::HighPrecisionStar
        );
        if high && !flagged && uncovered.is_some_and(|u| u <= 50.0) {
            out.set("status", "positive");
            out.set("message", format!("VNtyper calls {what}."));
        } else {
            let mut why = Vec::new();
            if flagged {
                why.push(Caveat::Flagged(flag.clone()));
            }
            if conf == Confidence::LowPrecision {
                why.push(Caveat::LowPrecision);
            }
            match uncovered {
                Some(u) if u > 50.0 => why.push(Caveat::Uncovered(u)),
                None => why.push(Caveat::NoUncoveredMeasure),
                Some(_) => {}
            }
            let why: Vec<String> = why.iter().map(Caveat::text).collect();
            out.set("status", "possible");
            out.set(
                "message",
                format!("VNtyper calls {what}, with caveats ({}).", why.join("; ")),
            );
        }
        return Ok(out);
    }
    if covered {
        out.set("status", "negative");
        out.set("message", "VNtyper calls no MUC1 VNTR frameshift.");
        out.set("advice", ADVICE_NONE);
        return Ok(out);
    }
    let depth = mean.map_or("no measured depth".into(), |m| {
        format!("{}x", crate::fixed(m, 0))
    });
    Ok(out.not_assessed(format!(
        "VNtyper calls no frameshift, but its coverage check is {qc} ({depth} over the VNTR, 50x is its floor on GRCh38), so a frameshift cannot be ruled out."
    )))
}

const SPEC: Spec = Spec {
    prog: "og vntyper report",
    usage: "usage: og vntyper report [-h] --run-dir RUN_DIR --out OUT [--exit-code EXIT_CODE] [--log LOG]\n",
    options: &["--run-dir", "--out", "--exit-code", "--log"],
    positionals: &[],
    required: &["--run-dir", "--out"],
};

/// `og vntyper VERB ARGS`: the exit code (2 for unreadable results or a
/// refused command line).
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if verb != Some("report") {
        return crate::unknown_verb("vntyper", "report", verb, err);
    }
    crate::run(&SPEC, rest, out, err, |a| {
        let exit_code = crate::int_opt(a, "--exit-code", 0)?;
        let run_dir = crate::opt(a, "--run-dir");
        Ok(report(
            run_dir,
            exit_code,
            a.get("--log"),
            crate::opt(a, "--out"),
        ))
    })
}

fn report(run_dir: &str, exit_code: i64, log: Option<&str>, out: &str) -> Outcome {
    let unreadable = |e: Error| {
        (
            Error(format!("cannot read VNtyper's results in {run_dir}: {e}")),
            2,
        )
    };
    let log_text = match log.filter(|l| Path::new(l).is_file()) {
        Some(l) => crate::read_lines(l, Utf8::Replace)
            .map_err(unreadable)?
            .concat(),
        None => String::new(),
    };
    let t = read(run_dir, exit_code, &log_text).map_err(unreadable)?;
    let mut kv = KvTable::metrics();
    for k in FIELDS {
        kv.push(k, t.get(k));
    }
    let mut body = Vec::new();
    kv.write_to(&mut body)
        .map_err(|e| (Error::io(&e, out), 2))?;
    crate::write_atomic(out, &body).map_err(|e| (e, 1))?;
    Ok(String::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    const KHEAD: &str = "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\tNomenclature\n";

    fn run(kestrel: Option<&str>, cov: Option<&str>, exit: i64, log: &str) -> Reading {
        let d = Dir::new("vntyper");
        if let Some(k) = kestrel {
            d.put("kestrel/kestrel_result.tsv", k);
        }
        if let Some(c) = cov {
            d.put("coverage/coverage_summary.tsv", c);
        }
        read(&d.root().display().to_string(), exit, log).unwrap()
    }

    fn cov(mean: &str, unc: &str, qc: &str) -> String {
        format!("mean\tpercent_uncovered\tcoverage_qc\n{mean}\t{unc}\t{qc}\n")
    }

    fn call(conf: &str, flag: &str) -> String {
        format!("## VNtyper\n{KHEAD}X\tdupC\t60\tG\tGC\tm\t12\t300\t0.04\t{conf}\t{flag}\tdupC\n")
    }

    #[test]
    fn statuses() {
        let r = run(
            Some(&call("High_Precision", "Not flagged")),
            Some(&cov("212.4", "0.0", "PASS")),
            0,
            "",
        );
        assert_eq!(r.get("status"), "positive");
        assert_eq!(r.get("vntr_mean_depth"), "212.40");
        let r = run(
            Some(&call("Low_Precision", "")),
            Some(&cov("212.4", "0.0", "PASS")),
            0,
            "",
        );
        assert_eq!(r.get("status"), "possible");
        assert!(r.get("message").contains("(low precision)"));
        let r = run(
            Some(&call("High_Precision*", "False_Positive_4bp")),
            Some(&cov("80", "60", "PASS")),
            0,
            "",
        );
        assert!(
            r.get("message").ends_with(
                "with caveats (flagged: False_Positive_4bp; 60% of the VNTR has no read)."
            )
        );
        let r = run(Some(KHEAD), Some(&cov("80", "0", "PASS")), 0, "");
        assert_eq!(r.get("status"), "not assessed");
        let r = run(
            Some(&format!(
                "{KHEAD}x\tNone\t\t\t\t\t\t\t\tNegative\tNot flagged\t\n"
            )),
            Some(&cov("30", "0", "FAIL")),
            0,
            "",
        );
        assert_eq!(r.get("status"), "not assessed");
        assert!(r.get("message").contains("is FAIL (30x over"));
        let r = run(None, None, 1, "x FASTQ layout 'empty' y");
        assert!(r.get("message").starts_with("The BAM has no read"));
        let r = run(None, Some(&cov("1", "1", "PASS")), 0, "");
        assert_eq!(
            r.get("message"),
            "VNtyper wrote no kestrel/kestrel_result.tsv."
        );
    }

    /// Red first: a call with exactly half the VNTR uncovered is positive.
    #[test]
    fn red_first_half_uncovered_is_still_positive() {
        let r = run(
            Some(&call("High_Precision", "Not flagged")),
            Some(&cov("212.4", "50", "PASS")),
            0,
            "",
        );
        assert_eq!(r.get("status"), "positive");
    }

    #[test]
    fn sensitivity_bands_and_ties() {
        assert!(
            sensitivity_line(Some(100.0))
                .unwrap()
                .contains("(at 100x or more)")
        );
        assert!(
            sensitivity_line(Some(50.0))
                .unwrap()
                .contains("(between 50x and 100x)")
        );
        // Half-even on the exact value, as Python's format.
        assert!(
            sensitivity_line(Some(48.5))
                .unwrap()
                .starts_with("Mean depth over the MUC1 VNTR: 48x.")
        );
        assert!(
            sensitivity_line(Some(49.5))
                .unwrap()
                .starts_with("Mean depth over the MUC1 VNTR: 50x.")
        );
        assert!(
            sensitivity_line(Some(10.0))
                .unwrap()
                .contains("in 48 of 109 runs (44%)")
        );
        assert!(sensitivity_line(Some(-1.0)).is_err());
    }
}
