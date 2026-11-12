//! `og y str-bed` and `og y str-report` (step 45): Y-STR repeat counts from
//! HipSTR's calls, in FamilyTreeDNA's marker order.
//!
//! ```text
//! og y str-bed --panel assets/ystr_panel.tsv --out regions.bed
//! og y str-report --panel assets/ystr_panel.tsv --vcf <id>_y_str.vcf.gz \
//!     --tsv <id>_y_str.tsv --csv <id>_y_str.csv [--log <id>_y_str.log] [--read-length 150]
//! ```
//!
//! `str-bed` writes HipSTR's region file: every region of the panel's
//! single and parts markers, and one region per span marker, sorted by
//! position. `str-report` writes one row per panel marker, in the panel's
//! order; with HipSTR's log, a region with no record carries the log's
//! reason.
//!
//! How a count is made: HipSTR writes per region GB, the base pairs the
//! call differs from the reference. A region's count is the panel's
//! reference copies (halves rounded up) plus GB / motif length, by floor
//! division, so GB -2 on 11 copies of a 4 bp motif is `10.2` (ten repeats
//! and two bases). Parts are summed; a span marker is one region; the
//! panel's offset is added last. Confidence: high (HipSTR's haploid filters
//! pass and 10 reads or more), low (the note says which filter failed), no
//! call, not callable.
//!
//! What changes on purpose: a motif length under 1, reference copies that
//! are negative or not finite, and an order that is not a whole number are
//! panel errors (the Python crashed later or divided by zero); an
//! unreadable VCF or log is an `ERROR:` line (exit 1), not a traceback.

use crate::Error;
use og_filters::args::Spec;
use og_io::Utf8;
use std::collections::HashMap;
use std::io::Write;

const MIN_READS: i128 = 10;
const MIN_Q: f64 = 0.9;
const MAX_STUTTER: f64 = 0.15;
const MAX_FLANK: f64 = 0.15;
const FLANK: i128 = 10;

pub const TSV_COLUMNS: [&str; 10] = [
    "order",
    "panel",
    "marker",
    "count",
    "confidence",
    "reads",
    "stutter_reads",
    "flank_indel_reads",
    "q",
    "note",
];
const NEED: [&str; 8] = [
    "order", "panel", "marker", "kind", "regions", "motif", "offset", "note",
];

/// One HipSTR region.
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    chrom: String,
    start: i128,
    end: i128,
    period: i128,
    copies: f64,
    name: String,
}

/// What a panel marker is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Single,
    Parts,
    Span,
    Multicopy,
    Absent,
}

/// One panel row.
#[derive(Debug, Clone)]
pub struct Marker {
    order: String,
    panel: String,
    marker: String,
    kind: Kind,
    /// The regions HipSTR genotypes: a span marker has one.
    regions: Vec<Region>,
    offset: Option<i128>,
    note: String,
}

/// The reference copies of a region, halves up.
fn halves_up(copies: f64) -> i128 {
    (copies + 0.5).trunc() as i128
}

/// The panel's rows, checked.
pub fn read_panel(path: &str) -> Result<Vec<Marker>, Error> {
    let lines: Vec<String> = crate::read_lines(path, Utf8::Strict)?
        .into_iter()
        .filter(|l| !crate::strip(l).is_empty() && !l.starts_with('#'))
        .map(|l| l.strip_suffix('\n').unwrap_or(&l).to_owned())
        .collect();
    let Some(head) = lines.first() else {
        return Err(Error(format!("{path}: no header")));
    };
    if head.split('\t').collect::<Vec<_>>() != NEED {
        return Err(Error(format!(
            "{path}: the header is not {}",
            NEED.join(" ")
        )));
    }
    let mut rows = Vec::new();
    for (i, line) in lines[1..].iter().enumerate() {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() != NEED.len() {
            return Err(Error(format!(
                "{path}: row {} has {} fields, not {}",
                i + 2,
                c.len(),
                NEED.len()
            )));
        }
        let marker = c[2];
        let kind = match c[3] {
            "single" => Kind::Single,
            "parts" => Kind::Parts,
            "span" => Kind::Span,
            "multicopy" => Kind::Multicopy,
            "absent" => Kind::Absent,
            k => return Err(Error(format!("{path}: {marker}: unknown kind {k}"))),
        };
        let mut regs = Vec::new();
        if c[4] != "." {
            for part in c[4].split(';') {
                let bad = || Error(format!("{path}: {marker}: unreadable region {part}"));
                let f: Vec<&str> = part.splitn(4, '/').collect();
                let [where_, period, copies, name] = f[..] else {
                    return Err(bad());
                };
                let w: Vec<&str> = where_.split(':').collect();
                let [chrom, span] = w[..] else {
                    return Err(bad());
                };
                let s: Vec<&str> = span.split('-').collect();
                let [start, end] = s[..] else {
                    return Err(bad());
                };
                let (Some(start), Some(end), Some(period), Some(copies)) = (
                    crate::int(start),
                    crate::int(end),
                    crate::int(period),
                    crate::float(copies),
                ) else {
                    return Err(bad());
                };
                if period < 1 || !copies.is_finite() || copies < 0.0 {
                    return Err(bad());
                }
                regs.push(Region {
                    chrom: chrom.into(),
                    start,
                    end,
                    period,
                    copies,
                    name: name.into(),
                });
            }
        }
        let has = matches!(kind, Kind::Single | Kind::Parts | Kind::Span);
        if has == regs.is_empty()
            || (kind == Kind::Single && regs.len() != 1)
            || (matches!(kind, Kind::Parts | Kind::Span) && regs.len() < 2)
        {
            return Err(Error(format!(
                "{path}: {marker}: kind {} with {} regions",
                c[3],
                regs.len()
            )));
        }
        if regs.iter().any(|r| r.period != regs[0].period) {
            return Err(Error(format!(
                "{path}: {marker}: parts with different motif lengths"
            )));
        }
        let offset = if c[6] == "." {
            None
        } else {
            Some(crate::int(c[6]).ok_or_else(|| {
                Error(format!(
                    "{path}: {marker}: offset {} is not a whole number",
                    c[6]
                ))
            })?)
        };
        if kind == Kind::Span {
            if regs.iter().any(|r| r.chrom != regs[0].chrom)
                || regs.windows(2).any(|w| w[1].start <= w[0].end)
            {
                return Err(Error(format!(
                    "{path}: {marker}: span parts not in order on one chromosome"
                )));
            }
            // One region over the whole locus, named after the marker; its
            // copies are the parts' (each rounded, halves up).
            let copies = regs.iter().map(|r| halves_up(r.copies)).sum::<i128>() as f64;
            regs = vec![Region {
                chrom: regs[0].chrom.clone(),
                start: regs[0].start,
                end: regs[regs.len() - 1].end,
                period: regs[0].period,
                copies,
                name: marker.into(),
            }];
        }
        rows.push(Marker {
            order: c[0].into(),
            panel: c[1].into(),
            marker: marker.into(),
            kind,
            regions: regs,
            offset,
            note: c[7].into(),
        });
    }
    for (i, r) in rows.iter().enumerate() {
        if crate::int(&r.order) != Some(i as i128 + 1) {
            return Err(Error(format!(
                "{path}: the order column is not 1 to {}",
                rows.len()
            )));
        }
    }
    Ok(rows)
}

/// HipSTR's region file.
pub fn bed(panel: &[Marker]) -> String {
    let mut seen: Vec<&Region> = Vec::new();
    let mut at: HashMap<&str, usize> = HashMap::new();
    for r in panel.iter().flat_map(|m| &m.regions) {
        match at.get(r.name.as_str()) {
            Some(&i) => seen[i] = r,
            None => {
                at.insert(&r.name, seen.len());
                seen.push(r);
            }
        }
    }
    seen.sort_by(|a, b| (&a.chrom, a.start, a.end).cmp(&(&b.chrom, b.start, b.end)));
    seen.iter()
        .map(|r| {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                r.chrom,
                r.start,
                r.end,
                r.period,
                pycompat::fmt::general(r.copies, 6),
                r.name
            )
        })
        .collect()
}

/// HipSTR's records by ID: the FORMAT values of the one sample.
fn read_calls(path: &str) -> Result<HashMap<String, HashMap<String, String>>, Error> {
    let mut calls = HashMap::new();
    let mut samples: Option<usize> = None;
    for line in crate::lines(path, Utf8::Replace)? {
        let line = line?;
        if line.starts_with("##") {
            continue;
        }
        let c: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(&line)
            .split('\t')
            .collect();
        if line.starts_with('#') {
            let n = c.len().saturating_sub(9);
            if n != 1 {
                return Err(Error(format!("{path}: {n} samples; step 45 reads one")));
            }
            samples = Some(n);
            continue;
        }
        if samples.is_none() || c.len() < 10 {
            return Err(Error(format!(
                "{path}: a record before the header, or one without a sample column"
            )));
        }
        let fmt = c[8]
            .split(':')
            .zip(c[9].split(':'))
            .map(|(k, v)| (k.to_owned(), v.to_owned()))
            .collect();
        calls.insert(c[2].to_owned(), fmt);
    }
    if samples.is_none() {
        return Err(Error(format!("{path}: no #CHROM header line")));
    }
    Ok(calls)
}

const LOG_REASONS: [&str; 5] = [
    "Skipping",
    "Aborting",
    "Stutter model training failed",
    "Haplotype construction failed",
    "Retraining stutter model training failed",
];

/// HipSTR's reason per region: (chrom, 1-based start, end) to the first
/// reason line after its "Processing region" line.
fn read_log(path: &str) -> Result<HashMap<(String, i128, i128), String>, Error> {
    let mut why = HashMap::new();
    let mut cur: Option<(String, i128, i128)> = None;
    for line in crate::read_lines(path, Utf8::Replace)? {
        let line = crate::strip(&line);
        if line.starts_with("Processing region ") {
            let c = pycompat::text::split(line, None, None);
            cur = match (
                c.get(2),
                c.get(3).and_then(|s| crate::int(s)),
                c.get(4).and_then(|s| crate::int(s)),
            ) {
                (Some(ch), Some(s), Some(e)) => Some(((*ch).to_owned(), s + 1, e)),
                _ => None,
            };
        } else if let Some(k) = &cur
            && !why.contains_key(k)
            && LOG_REASONS.iter().any(|r| line.starts_with(r))
        {
            why.insert(k.clone(), line.to_owned());
        }
    }
    Ok(why)
}

/// A called region: bp difference, reads, stutter, flank indels, Q, and the
/// filters it fails.
struct PartCall {
    bp: i128,
    dp: i128,
    st: i128,
    fl: i128,
    q: Option<f64>,
    problems: Vec<String>,
}

fn part_call(fmt: &HashMap<String, String>) -> Result<PartCall, String> {
    let get = |k: &str| fmt.get(k).map(String::as_str);
    let (gt, gb) = (get("GT").unwrap_or("."), get("GB").unwrap_or("."));
    let none = |s: &str| matches!(s, "." | "./." | ".|.");
    if none(gt) || none(gb) {
        let why = get("FILTER").filter(|w| !w.is_empty() && *w != "." && *w != "PASS");
        return Err(format!(
            "HipSTR made no genotype{}",
            why.map_or(String::new(), |w| format!(" ({w})"))
        ));
    }
    let mut lens: Vec<Option<i128>> = gb.replace('|', "/").split('/').map(crate::int).collect();
    lens.sort();
    lens.dedup();
    if lens.contains(&None) {
        return Err(format!("unreadable GB {gb}"));
    }
    if lens.len() > 1 {
        return Err(format!("two lengths on the Y chromosome (GB {gb})"));
    }
    let bp = lens[0].unwrap_or(0);
    let n = |k: &str| get(k).and_then(crate::int).unwrap_or(0);
    let (dp, st, fl) = (n("DP"), n("DSTUTTER"), n("DFLANKINDEL"));
    let q = get("Q").and_then(crate::float);
    let mut problems = Vec::new();
    if q.is_none_or(|q| q < MIN_Q) {
        problems.push(format!(
            "posterior {} below {}",
            get("Q").unwrap_or("."),
            pycompat::fmt::repr(MIN_Q)
        ));
    }
    if dp != 0 && st as f64 / dp as f64 > MAX_STUTTER {
        problems.push(format!("{st} of {dp} reads with stutter"));
    }
    if dp != 0 && fl as f64 / dp as f64 > MAX_FLANK {
        problems.push(format!("{fl} of {dp} reads with a flank indel"));
    }
    if dp < MIN_READS {
        problems.push(format!("{dp} reads, fewer than {MIN_READS}"));
    }
    Ok(PartCall {
        bp,
        dp,
        st,
        fl,
        q,
        problems,
    })
}

/// Python's `divmod` on ints: the quotient rounded down.
fn divmod(a: i128, b: i128) -> (i128, i128) {
    let (mut q, mut r) = (a / b, a % b);
    if r != 0 && ((r < 0) != (b < 0)) {
        q -= 1;
        r += b;
    }
    (q, r)
}

/// The count: reference copies plus the bp difference in repeats, the rest
/// after a dot.
fn count_text(ref_copies: i128, bp: i128, period: i128, offset: Option<i128>) -> String {
    let (whole, rest) = divmod(bp, period);
    let n = ref_copies + whole + offset.unwrap_or(0);
    if rest != 0 {
        format!("{n}.{rest}")
    } else {
        n.to_string()
    }
}

/// `min` as Python's builtin: the first item no later one is smaller than.
fn py_min(xs: &[f64]) -> f64 {
    xs[1..]
        .iter()
        .fold(xs[0], |m, &x| if x < m { x } else { m })
}

/// One marker's row, by [`TSV_COLUMNS`].
pub fn marker_row(
    m: &Marker,
    calls: &HashMap<String, HashMap<String, String>>,
    read_length: i64,
    log: Option<&HashMap<(String, i128, i128), String>>,
) -> Vec<String> {
    let mut row: Vec<String> = vec![
        m.order.clone(),
        m.panel.clone(),
        m.marker.clone(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
    ];
    let notes: Vec<String> = if m.note == "." {
        vec![]
    } else {
        vec![m.note.clone()]
    };
    let finish = |mut row: Vec<String>, conf: &str, note: String| {
        row[4] = conf.into();
        row[9] = note;
        row
    };
    match m.kind {
        Kind::Multicopy => {
            return finish(
                row,
                "not callable",
                "not callable from short reads: a multi-copy marker in duplicated Y sequence"
                    .into(),
            );
        }
        Kind::Absent => {
            return finish(
                row,
                "not callable",
                "not callable: no region for it in the HipSTR marker set".into(),
            );
        }
        _ => {}
    }
    let longest = m
        .regions
        .iter()
        .map(|r| r.end - r.start + 1)
        .max()
        .unwrap_or(0);
    if longest + 2 * FLANK > i128::from(read_length) {
        return finish(
            row,
            "not callable",
            format!(
                "not callable from short reads: a {longest} bp repeat is longer than {read_length} bp reads can span with {FLANK} bases on each side"
            ),
        );
    }
    let period = m.regions[0].period;
    let (mut bp_sum, mut ref_sum) = (0i128, 0i128);
    let (mut dps, mut sts, mut fls, mut qs) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut problems = Vec::new();
    let many = m.regions.len() > 1;
    for r in &m.regions {
        let Some(fmt) = calls.get(&r.name) else {
            let why = match log {
                None => format!("HipSTR wrote no record for {}; see HipSTR's log", r.name),
                Some(l) => match l.get(&(r.chrom.clone(), r.start, r.end)) {
                    Some(reason) => {
                        format!("HipSTR wrote no record for {}; its log: {reason}", r.name)
                    }
                    None => format!(
                        "HipSTR wrote no record for {}, and its log gives no reason",
                        r.name
                    ),
                },
            };
            let mut n = notes.clone();
            n.push(why);
            return finish(row, "no call", n.join("; "));
        };
        let got = match part_call(fmt) {
            Ok(g) => g,
            Err(why) => {
                let mut n = notes.clone();
                n.push(format!("{}: {why}", r.name));
                return finish(row, "no call", n.join("; "));
            }
        };
        bp_sum += got.bp;
        ref_sum += halves_up(r.copies);
        dps.push(got.dp);
        sts.push(got.st);
        fls.push(got.fl);
        qs.push(got.q);
        problems.extend(
            got.problems
                .into_iter()
                .map(|x| if many { format!("{}: {x}", r.name) } else { x }),
        );
    }
    row[3] = count_text(ref_sum, bp_sum, period, m.offset);
    row[4] = if problems.is_empty() { "high" } else { "low" }.into();
    row[5] = dps.iter().min().map_or(String::new(), ToString::to_string);
    row[6] = sts.iter().max().map_or(String::new(), ToString::to_string);
    row[7] = fls.iter().max().map_or(String::new(), ToString::to_string);
    row[8] = match qs.iter().copied().collect::<Option<Vec<f64>>>() {
        Some(q) => pycompat::fmt::general(py_min(&q), 6),
        None => ".".into(),
    };
    let mut n = notes;
    n.extend(problems);
    row[9] = if n.is_empty() {
        ".".into()
    } else {
        n.join("; ")
    };
    row
}

/// The files of `og y str-report`.
pub struct Files<'a> {
    pub panel: &'a str,
    pub vcf: &'a str,
    pub tsv: &'a str,
    pub csv: &'a str,
    pub log: Option<&'a str>,
}

/// `og y str-report`: writes both tables; the stdout line.
pub fn report(f: &Files<'_>, read_length: i64) -> Result<String, Error> {
    let panel = read_panel(f.panel)?;
    let calls = read_calls(f.vcf)?;
    let log = match f.log {
        Some(l) => Some(read_log(l)?),
        None => None,
    };
    let rows: Vec<Vec<String>> = panel
        .iter()
        .map(|m| marker_row(m, &calls, read_length, log.as_ref()))
        .collect();
    let mut tsv = TSV_COLUMNS.join("\t") + "\n";
    for r in &rows {
        tsv.push_str(
            &r.iter()
                .map(|c| c.replace('\t', " "))
                .collect::<Vec<_>>()
                .join("\t"),
        );
        tsv.push('\n');
    }
    let mut w = pycompat::csv::Writer::new(Vec::new(), ',');
    let io = |e: std::io::Error| Error::io(&e, f.csv);
    w.write_row(["marker", "value", "confidence"]).map_err(io)?;
    for r in &rows {
        w.write_row([r[2].as_str(), &r[3], &r[4]]).map_err(io)?;
    }
    crate::write_atomic(f.tsv, tsv.as_bytes())?;
    crate::write_atomic(f.csv, &w.into_inner())?;
    let tally = |k: &str| rows.iter().filter(|r| r[4] == k).count();
    Ok(format!(
        "Y-STR markers: {} high, {} low, {} no call, {} not callable\n",
        tally("high"),
        tally("low"),
        tally("no call"),
        tally("not callable")
    ))
}

const BED: Spec = Spec {
    prog: "og y str-bed",
    usage: "usage: og y str-bed [-h] --panel PANEL --out OUT\n",
    options: &["--panel", "--out"],
    positionals: &[],
    required: &["--panel", "--out"],
};

const REPORT: Spec = Spec {
    prog: "og y str-report",
    usage: "usage: og y str-report [-h] --panel PANEL --vcf VCF --tsv TSV --csv CSV [--log LOG] [--read-length READ_LENGTH]\n",
    options: &[
        "--panel",
        "--vcf",
        "--tsv",
        "--csv",
        "--log",
        "--read-length",
    ],
    positionals: &[],
    required: &["--panel", "--vcf", "--tsv", "--csv"],
};

/// `og y VERB ARGS`: the exit code (1 for a problem, 2 for a refused
/// command line).
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    match verb {
        Some("str-bed") => crate::run(&BED, rest, out, err, |a| {
            Ok((|| {
                let panel = read_panel(crate::opt(a, "--panel"))?;
                crate::write_atomic(crate::opt(a, "--out"), bed(&panel).as_bytes())?;
                Ok(String::new())
            })()
            .map_err(|e: Error| (e, 1)))
        }),
        Some("str-report") => crate::run(&REPORT, rest, out, err, |a| {
            let read_length = crate::int_opt(a, "--read-length", 150)?;
            let files = Files {
                panel: crate::opt(a, "--panel"),
                vcf: crate::opt(a, "--vcf"),
                tsv: crate::opt(a, "--tsv"),
                csv: crate::opt(a, "--csv"),
                log: a.get("--log"),
            };
            Ok(report(&files, read_length).map_err(|e| (e, 1)))
        }),
        other => crate::unknown_verb("y", "str-bed|str-report", other, err),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floor_division_of_a_negative_difference() {
        assert_eq!(count_text(11, -2, 4, None), "10.2");
        assert_eq!(count_text(11, 0, 4, None), "11");
        assert_eq!(count_text(11, 6, 4, Some(-1)), "11.2");
        assert_eq!(divmod(-7, 4), (-2, 1));
        assert_eq!(divmod(7, -2), (-4, -1));
    }

    /// Red first: GB -4 on a 4 bp motif is one repeat fewer, not `11.-4`.
    #[test]
    fn red_first_truncating_division_would_be_wrong() {
        assert_eq!(count_text(11, -4, 4, None), "10");
        assert_eq!(count_text(11, -1, 4, None), "10.3");
    }

    #[test]
    fn part_calls() {
        let f = |s: &str| -> HashMap<String, String> {
            s.split(';')
                .map(|kv| {
                    let (k, v) = kv.split_once('=').unwrap();
                    (k.into(), v.into())
                })
                .collect()
        };
        let p = part_call(&f("GT=1;GB=-4;Q=0.95;DP=30;DSTUTTER=6;DFLANKINDEL=0")).unwrap();
        assert_eq!(p.problems, ["6 of 30 reads with stutter"]);
        let p = part_call(&f("GT=1;GB=0;Q=0.5;DP=5;DSTUTTER=0;DFLANKINDEL=0")).unwrap();
        assert_eq!(
            p.problems,
            ["posterior 0.5 below 0.9", "5 reads, fewer than 10"]
        );
        assert_eq!(
            part_call(&f("GT=.;GB=.;FILTER=LowQ")).err().unwrap(),
            "HipSTR made no genotype (LowQ)"
        );
        assert_eq!(
            part_call(&f("GT=0/1;GB=0/4")).err().unwrap(),
            "two lengths on the Y chromosome (GB 0/4)"
        );
        assert_eq!(part_call(&f("GT=0;GB=x")).err().unwrap(), "unreadable GB x");
        assert!(py_min(&[f64::NAN, 0.5]).is_nan());
        assert_eq!(py_min(&[0.5, f64::NAN]), 0.5);
    }
}
