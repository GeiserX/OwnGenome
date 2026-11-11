//! `og cn <verb>`:
//!
//! ```text
//! og cn bed [--catalog TSV] [--locus L ...] [--kind interval] [--psv]
//! og cn call [--catalog TSV] --locus L --bedcov T [--bedcov-mapq1 T1] [--pileup P]
//!            [--sex male|female|unknown] [--table OUT.tsv] [--bam-header H]
//!            [--marker LPA_KIV2_4925 --marker-pileup P --marker-out OUT.tsv]
//! og cn intervals [--catalog TSV] --bedcov T --sex male|female|unknown --out OUT.tsv
//! og cn intervals [--catalog TSV] --failed --out OUT.tsv
//! og cn check-bam --header H
//! og cn marker-regions NAME
//! og cn marker-read NAME --pileup P [--copy-table T] --out OUT.tsv [--catalog TSV]
//! ```
//!
//! The commands of `copy_number.py` with the same options (`--json`, which
//! no step passed, is gone), `traits_report.py --marker-regions` and
//! `--marker`, and two additions: `call --sex`, so a gene row on chrX or
//! chrY is read against the sample's count there (the Python never passed a
//! sex, so such a row could only read `sex_unknown`), and `call --marker`,
//! which reads the KIV-2 4925G>A marker beside the count in the same run and
//! judges the aligner on the count's own MAPQ 0 share. `marker-read
//! LPA_KIV2_4925 --copy-table` reads a table written earlier, where the
//! share is printed at 3 decimals, and judges that printed share.
//!
//! Exit codes: `call` 0 whatever the status, 2 on unreadable input, 3 on a
//! broken control, 4 on a BAM with ALT contigs; `check-bam` 0 or 4 (2 when
//! the header cannot be read); `bed`, `marker-regions` and `intervals` 0
//! or 2 (2 also when stdout or the table cannot be written); `marker-read`
//! 0 whatever the reading (an unreadable pileup reads 'not assessed'), 2
//! when the table cannot be written. Usage errors exit 2.

use crate::catalog::{self, Catalog, Input, Locus};
use crate::engine::{self, AssessError, Block, Means, Reason, Sex, Status};
use crate::layout::{self, INTERVAL_COLUMNS};
use crate::markers::{self, LPA_KIV2_4925, Marker, SMN1_27134};
use crate::{Error, input, smn};
use og_filters::py_repr;
use std::collections::HashMap;
use std::io::Write;

const USAGE: &str = "usage: og cn bed [--catalog TSV] [--locus L ...] [--kind interval] [--psv]
       og cn call [--catalog TSV] --locus L --bedcov T [--bedcov-mapq1 T1] [--pileup P]
                  [--sex male|female|unknown] [--table OUT] [--bam-header H]
                  [--marker LPA_KIV2_4925 --marker-pileup P --marker-out OUT]
       og cn intervals [--catalog TSV] (--bedcov T --sex male|female|unknown | --failed) --out OUT
       og cn check-bam --header H
       og cn marker-regions NAME
       og cn marker-read NAME --pileup P [--copy-table T] --out OUT [--catalog TSV]
";

/// A parsed command line: options with a value (`--locus` may repeat),
/// flags, positionals.
#[derive(Debug, Default)]
struct Args {
    opts: HashMap<&'static str, Vec<String>>,
    flags: Vec<&'static str>,
    pos: Vec<String>,
}

impl Args {
    fn get(&self, name: &str) -> Option<&str> {
        self.opts
            .get(name)
            .and_then(|v| v.last())
            .map(String::as_str)
    }

    fn all(&self, name: &str) -> &[String] {
        self.opts.get(name).map_or(&[], Vec::as_slice)
    }

    fn flag(&self, name: &str) -> bool {
        self.flags.contains(&name)
    }
}

/// Parses `args` by the options a verb takes, as argparse does: `--name
/// VALUE` or `--name=VALUE`, flags without a value, positionals in order.
fn parse(
    args: &[String],
    options: &[&'static str],
    flags: &[&'static str],
    positionals: usize,
    required: &[&str],
) -> Result<Args, String> {
    let mut out = Args::default();
    let mut extra = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        i += 1;
        if a == "-h" || a == "--help" {
            return Err(String::new());
        }
        if a == "-" || !a.starts_with('-') {
            if out.pos.len() < positionals {
                out.pos.push(a.clone());
            } else {
                extra.push(a.clone());
            }
            continue;
        }
        let (name, inline) = match a.split_once('=') {
            Some((n, v)) => (n, Some(v.to_owned())),
            None => (a.as_str(), None),
        };
        if let Some(&f) = flags.iter().find(|f| **f == name) {
            if inline.is_some() {
                return Err(format!(
                    "argument {f}: ignored explicit argument {}",
                    py_repr(inline.as_deref().unwrap_or(""))
                ));
            }
            out.flags.push(f);
            continue;
        }
        let Some(&known) = options.iter().find(|o| **o == name) else {
            extra.push(a.clone());
            continue;
        };
        let value = match inline {
            Some(v) => v,
            None => match args.get(i) {
                Some(v) if v == "-" || !v.starts_with('-') => {
                    i += 1;
                    v.clone()
                }
                _ => return Err(format!("argument {known}: expected one argument")),
            },
        };
        out.opts.entry(known).or_default().push(value);
    }
    if out.pos.len() < positionals {
        return Err("the following arguments are required: NAME".into());
    }
    let missing: Vec<&str> = required
        .iter()
        .copied()
        .filter(|r| !out.opts.contains_key(r))
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    if !extra.is_empty() {
        return Err(format!("unrecognized arguments: {}", extra.join(" ")));
    }
    Ok(out)
}

fn usage(verb: &str, msg: &str, err: &mut dyn Write) -> u8 {
    let _ = write!(err, "{USAGE}");
    let _ = writeln!(err, "og cn {verb}: error: {msg}");
    2
}

/// `og cn <verb> [args]`: the exit code.
pub fn main(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let Some(verb) = args.first() else {
        return usage("", "the following arguments are required: verb", err);
    };
    let rest = &args[1..];
    let spec: (&[&'static str], &[&'static str], usize, &[&str]) = match verb.as_str() {
        "bed" => (&["--catalog", "--locus", "--kind"], &["--psv"], 0, &[]),
        "call" => (
            &[
                "--catalog",
                "--locus",
                "--bedcov",
                "--bedcov-mapq1",
                "--table",
                "--bam-header",
                "--pileup",
                "--sex",
                "--marker",
                "--marker-pileup",
                "--marker-out",
            ],
            &[],
            0,
            &["--locus", "--bedcov"],
        ),
        "intervals" => (
            &["--catalog", "--bedcov", "--sex", "--out"],
            &["--failed"],
            0,
            &["--out"],
        ),
        "check-bam" => (&["--header"], &[], 0, &["--header"]),
        "marker-regions" => (&[], &[], 1, &[]),
        "marker-read" => (
            &["--pileup", "--copy-table", "--out", "--catalog"],
            &[],
            1,
            &["--pileup", "--out"],
        ),
        "-h" | "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        other => return usage("", &format!("unknown cn verb {}", py_repr(other)), err),
    };
    let a = match parse(rest, spec.0, spec.1, spec.2, spec.3) {
        Ok(a) => a,
        Err(m) if m.is_empty() => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        Err(m) => return usage(verb, &m, err),
    };
    if let Some(m) = choice_error(&a) {
        return usage(verb, &m, err);
    }
    let code = match verb.as_str() {
        "bed" => bed(&a, out),
        "call" => call(&a, out, err),
        "intervals" => intervals(&a, out, err),
        "check-bam" => check_bam(&a, err),
        "marker-regions" => {
            let m = markers::marker(&a.pos[0]).expect("checked");
            markers::regions(m)
                .iter()
                .try_for_each(|r| writeln!(out, "{r}"))
                .and_then(|()| out.flush())
                .map_err(unwritten)
                .map(|()| 0)
        }
        _ => marker_read(&a, out, err),
    };
    let _ = out.flush();
    match code {
        Ok(c) => c,
        Err((c, e)) => {
            let _ = writeln!(err, "ERROR: {e}");
            c
        }
    }
}

/// The argparse `choices` checks.
fn choice_error(a: &Args) -> Option<String> {
    let bad = |name: &str, v: &str, choices: &[&str]| {
        Some(format!(
            "argument {name}: invalid choice: {} (choose from {})",
            py_repr(v),
            choices
                .iter()
                .map(|c| py_repr(c))
                .collect::<Vec<_>>()
                .join(", ")
        ))
    };
    if let Some(k) = a.get("--kind").filter(|k| *k != "interval") {
        return bad("--kind", k, &["interval"]);
    }
    if let Some(s) = a.get("--sex").filter(|s| Sex::parse(s).is_none()) {
        return bad("--sex", s, &["male", "female", "unknown"]);
    }
    if let Some(m) = a.get("--marker").filter(|m| *m != LPA_KIV2_4925.key) {
        return bad("--marker", m, &[LPA_KIV2_4925.key]);
    }
    if let Some(n) = a.pos.first().filter(|n| markers::marker(n).is_none()) {
        return bad("NAME", n, &[LPA_KIV2_4925.key, SMN1_27134.key]);
    }
    None
}

type Fail = (u8, Error);

fn fail2(e: Error) -> Fail {
    (2, e)
}

/// A region list that could not be written in full: exit 2, so a step
/// never reads a truncated list.
fn unwritten(e: std::io::Error) -> Fail {
    fail2(Error(format!("cannot write the regions: {e}")))
}

fn load(a: &Args) -> Result<Catalog, Error> {
    match a.get("--catalog") {
        Some(p) => Catalog::load(p),
        None => {
            let p = catalog::default_path()
                .ok_or_else(|| Error("no bin/cn_loci.tsv beside og; give --catalog".to_owned()))?;
            Catalog::load(&p.to_string_lossy())
        }
    }
}

fn bed(a: &Args, out: &mut dyn Write) -> Result<u8, Fail> {
    let cat = load(a).map_err(fail2)?;
    let loci = a.all("--locus");
    if a.flag("--psv") {
        if loci.is_empty() || a.get("--kind").is_some() {
            return Err(fail2(Error(
                "bed --psv needs --locus (and no --kind)".into(),
            )));
        }
        for loc in loci {
            let Some(r) = cat.loci.get(loc).filter(|r| !r.psv.is_empty()) else {
                return Err(fail2(Error(format!(
                    "no psv locus {} in the catalog",
                    py_repr(loc)
                ))));
            };
            let mut pos: Vec<(i64, String)> = r
                .psv
                .iter()
                .map(|p| (p.gene_pos, p.name.clone()))
                .chain(
                    r.psv
                        .iter()
                        .map(|p| (p.pseudo_pos, format!("{}_pseudo", p.name))),
                )
                .collect();
            pos.sort();
            for (s, n) in pos {
                writeln!(out, "{}\t{s}\t{}\t{n}", r.chrom, s + 1).map_err(unwritten)?;
            }
        }
        out.flush().map_err(unwritten)?;
        return Ok(0);
    }
    let names: Vec<String> = if a.get("--kind").is_some() {
        if !loci.is_empty() {
            return Err(fail2(Error(
                "bed: --kind interval or --locus, not both".into(),
            )));
        }
        cat.loci
            .values()
            .filter(|r| r.input == Input::Bedcov && r.is_interval_row())
            .map(|r| r.id.clone())
            .collect()
    } else if loci.is_empty() {
        cat.loci
            .values()
            .filter(|r| r.input == Input::Bedcov)
            .map(|r| r.id.clone())
            .collect()
    } else {
        loci.to_vec()
    };
    for loc in &names {
        let Some(r) = cat.bedcov(loc) else {
            return Err(fail2(Error(format!(
                "no bedcov locus {} in the catalog",
                py_repr(loc)
            ))));
        };
        for (c, s, e, n) in r.regions() {
            writeln!(out, "{c}\t{s}\t{e}\t{n}").map_err(unwritten)?;
        }
    }
    out.flush().map_err(unwritten)?;
    Ok(0)
}

fn check_bam(a: &Args, err: &mut dyn Write) -> Result<u8, Fail> {
    match input::check_bam(a.get("--header").expect("required")).map_err(fail2)? {
        Some(refusal) => {
            let _ = writeln!(err, "{refusal}");
            Ok(4)
        }
        None => Ok(0),
    }
}

/// The count of one locus, and with `--marker`, the KIV-2 4925G>A reading.
fn call(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> Result<u8, Fail> {
    if let Some(h) = a.get("--bam-header")
        && let Some(refusal) = input::check_bam(h).map_err(fail2)?
    {
        let _ = writeln!(err, "{refusal}");
        return Ok(4);
    }
    let marker = match (
        a.get("--marker"),
        a.get("--marker-pileup"),
        a.get("--marker-out"),
    ) {
        (None, None, None) => None,
        (Some(_), Some(p), Some(o)) if a.get("--locus") == Some("LPA_KIV2") => Some((p, o)),
        _ => {
            return Err(fail2(Error(
                "--marker LPA_KIV2_4925 goes with --locus LPA_KIV2, --marker-pileup and --marker-out".into(),
            )));
        }
    };
    let sex = a.get("--sex").and_then(Sex::parse).unwrap_or_default();
    let cat = load(a).map_err(fail2)?;
    let locus = a.get("--locus").expect("required");
    let row = cat.bedcov(locus).ok_or_else(|| {
        fail2(Error(format!(
            "no bedcov locus {} in the catalog",
            py_repr(locus)
        )))
    })?;
    let block = count(a, row, sex).map_err(|e| match e {
        Count::Broken(s) => (3, Error(s)),
        Count::Input(e) => (2, e),
    })?;
    let (rows, line) = layout::table(row, &block);
    if let Some(t) = a.get("--table") {
        write_table(t, &rows).map_err(|e| (2, Error(format!("cannot write: {e}"))))?;
    }
    let _ = writeln!(out, "{line}");
    if let Some((pileup, path)) = marker {
        let kiv2 = rows
            .iter()
            .find(|(k, _)| k == "kiv2_copies")
            .map_or(".", |(_, v)| v.as_str());
        let why = match (block.mapq0_share, row.depth().mapq0_min) {
            (Some(s), Some(min)) => markers::mapq0_why_not(s, min),
            _ => None,
        };
        let counts = marker_counts(pileup, &LPA_KIV2_4925, "KIV-2 4925G>A is not assessed", err);
        let rows = markers::null_allele(counts, kiv2, why.as_deref());
        write_table(path, &rows).map_err(|e| (2, Error(format!("cannot write: {e}"))))?;
        let _ = writeln!(out, "LPA {}", rows[8].1);
    }
    Ok(0)
}

enum Count {
    Broken(String),
    Input(Error),
}

impl From<Error> for Count {
    fn from(e: Error) -> Count {
        Count::Input(e)
    }
}

fn count(a: &Args, row: &Locus, sex: Sex) -> Result<Block, Count> {
    let bedcov = a.get("--bedcov").expect("required");
    let m = input::read_bedcov(bedcov, row)?;
    let dp = row.depth();
    let q1_path = a.get("--bedcov-mapq1");
    if (dp.mapq0_min.is_some() || dp.mapq0_max.is_some()) && q1_path.is_none() {
        return Err(Error(format!(
            "{} has a MAPQ 0 guard: --bedcov-mapq1 is needed",
            row.id
        ))
        .into());
    }
    let q1 = q1_path.map(|p| input::read_bedcov(p, row)).transpose()?;
    if let (Some(q), Some(p)) = (&q1, q1_path)
        && m.is_empty() != q.is_empty()
    {
        let which = if m.is_empty() { bedcov } else { p };
        return Err(Error(format!(
                "{which} is empty but the other table is not: both must cover the regions, or both be empty (no reads on that contig)"
            ))
            .into());
    }
    let block = engine::assess(row, &m, q1.as_ref(), sex).map_err(|e| match e {
        AssessError::BrokenControl(s) => Count::Broken(s),
        AssessError::NeedsMapq1(s) => Count::Input(Error(s)),
    })?;
    if row.psv.is_empty() {
        return Ok(block);
    }
    let Some(pileup) = a.get("--pileup") else {
        return Err(Error(format!(
            "{} is a psv row: --pileup (samtools mpileup over `og cn bed --psv`) is needed",
            row.id
        ))
        .into());
    };
    let counts = if m.is_empty() {
        HashMap::new()
    } else {
        input::read_psv_pileup(pileup, row)?
    };
    Ok(engine::psv_assess(row, block, &counts))
}

fn intervals(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> Result<u8, Fail> {
    let cat = load(a).map_err(fail2)?;
    let rows: Vec<&Locus> = cat.loci.values().filter(|r| r.is_interval_row()).collect();
    let lines: Vec<Vec<String>> = if a.flag("--failed") {
        rows.iter()
            .map(|r| {
                let mut b = Block::empty();
                b.reason = Some(Reason::CheckFailed);
                layout::interval_line(r, &b)
            })
            .collect()
    } else {
        let (Some(bedcov), Some(sex)) = (a.get("--bedcov"), a.get("--sex").and_then(Sex::parse))
        else {
            return Err(fail2(Error(
                "intervals needs --bedcov and --sex (or --failed)".into(),
            )));
        };
        call_intervals(
            &rows,
            &input::read_bedcov_all(bedcov).map_err(fail2)?,
            sex,
            err,
        )
        .map_err(fail2)?
    };
    let path = a.get("--out").expect("required");
    let mut body = INTERVAL_COLUMNS.join("\t");
    body.push('\n');
    for l in &lines {
        body.push_str(&l.join("\t"));
        body.push('\n');
    }
    write_atomic(path, body.as_bytes()).map_err(|e| (2, Error(format!("cannot write: {e}"))))?;
    let mut n: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
    for l in &lines {
        *n.entry(l[5].as_str()).or_default() += 1;
    }
    let parts: Vec<String> = n.iter().map(|(k, v)| format!("{v} {k}")).collect();
    let _ = writeln!(out, "interval rows: {} ({})", lines.len(), parts.join(", "));
    Ok(0)
}

/// Every interval row from one bedcov table; a broken control makes that
/// row `not assessed` (broken_control), never the whole table.
fn call_intervals(
    rows: &[&Locus],
    table: &[input::BedcovLine],
    sex: Sex,
    err: &mut dyn Write,
) -> Result<Vec<Vec<String>>, Error> {
    let mut by_name: HashMap<&str, &input::BedcovLine> = HashMap::new();
    for l in table {
        if by_name.insert(&l.name, l).is_some() {
            return Err(Error(format!(
                "segment {} appears twice in the bedcov table",
                l.name
            )));
        }
    }
    let mut out = Vec::new();
    for row in rows {
        let want = row.regions();
        if want.iter().any(|w| !by_name.contains_key(w.3)) {
            return Err(Error(format!(
                "the bedcov table lacks some of {}'s segments: its regions differ from `og cn bed --kind interval`",
                row.id
            )));
        }
        let mut m = Means::new();
        for (c, s, e, n) in want {
            let l = by_name[n];
            if (l.chrom.as_str(), l.start, l.end) != (c, s, e) {
                return Err(Error(format!(
                    "segment {n} of {} has other coordinates in the bedcov table",
                    row.id
                )));
            }
            m.insert(n.to_owned(), l.sum / (e - s) as f64);
        }
        let block = match engine::assess(row, &m, None, sex) {
            Ok(b) => b,
            Err(e) => {
                let _ = writeln!(err, "WARNING: {e}");
                let mut b = engine::assess(row, &Means::new(), None, sex)
                    .map_err(|e| Error(e.to_string()))?;
                if b.status == Status::NotAssessed {
                    b.reason = Some(Reason::BrokenControl);
                }
                b
            }
        };
        out.push(layout::interval_line(row, &block));
    }
    Ok(out)
}

/// Reads a marker's pileup; an unreadable one is a warning and no counts.
fn marker_counts(
    pileup: &str,
    m: &Marker,
    what: &str,
    err: &mut dyn Write,
) -> Option<markers::Counts> {
    match markers::read_pileup(pileup, m) {
        Ok(c) => Some(c),
        Err(e) => {
            let _ = writeln!(err, "WARNING: {e}; {what}");
            None
        }
    }
}

fn marker_read(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> Result<u8, Fail> {
    let m = markers::marker(&a.pos[0]).expect("checked");
    let pileup = a.get("--pileup").expect("required");
    let path = a.get("--out").expect("required");
    let copy_table = a.get("--copy-table");
    let (rows, line) = if LPA_KIV2_4925.key == m.key {
        // Step 46's table written earlier: the count, and the MAPQ 0 share
        // as it was printed (3 decimals).
        let t: HashMap<String, String> = crate::read_lines(copy_table.unwrap_or(""))
            .map(|ls| {
                ls.iter()
                    .filter_map(|l| l.split_once('\t'))
                    .map(|(k, v)| (k.to_owned(), v.to_owned()))
                    .collect()
            })
            .unwrap_or_default();
        let kiv2 = t.get("kiv2_copies").map_or(".", String::as_str);
        let why = match t.get("mapq0_share").map(|s| pycompat::num::parse_float(s)) {
            Some(Ok(share)) => {
                let cat = load(a).map_err(fail2)?;
                let floor = cat
                    .bedcov("LPA_KIV2")
                    .and_then(|r| r.depth().mapq0_min)
                    .ok_or_else(|| {
                        fail2(Error(
                            "no LPA_KIV2 row with a MAPQ 0 floor in the catalog".into(),
                        ))
                    })?;
                markers::mapq0_why_not(share, floor)
            }
            _ => None,
        };
        let counts = marker_counts(pileup, m, "KIV-2 4925G>A is not assessed", err);
        let rows = markers::null_allele(counts, kiv2, why.as_deref());
        let line = format!("LPA {}", rows[8].1);
        (rows, line)
    } else {
        let counts = marker_counts(pileup, m, "the marker is not assessed", err);
        let read = match smn::read(copy_table, m.chrom, m.copies) {
            Ok(r) => r,
            Err(e) => {
                let _ = writeln!(err, "WARNING: {e}; the SMN copy numbers are not read");
                smn::SmnRead {
                    pscn: None,
                    text: "not read (the SMN table could not be read)".into(),
                    agcn: None,
                }
            }
        };
        let rows = markers::smn1_marker(counts, &read);
        let line = rows[9].1.clone();
        (rows, line)
    };
    write_table(path, &rows).map_err(|e| (2, Error(format!("cannot write: {e}"))))?;
    let _ = writeln!(out, "{line}");
    Ok(0)
}

/// Writes a `metric<TAB>value` table whole, through a temporary file.
fn write_table(path: &str, rows: &[(String, String)]) -> Result<(), String> {
    let mut body = String::from("metric\tvalue\n");
    for (k, v) in rows {
        body.push_str(k);
        body.push('\t');
        body.push_str(v);
        body.push('\n');
    }
    write_atomic(path, body.as_bytes())
}

fn write_atomic(path: &str, body: &[u8]) -> Result<(), String> {
    let run = || -> std::io::Result<()> {
        let mut f = og_io::AtomicFile::create(path)?;
        f.write_all(body)?;
        f.commit()
    };
    run().map_err(|e| og_filters::os_error(&e, path))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let c = main(&args, &mut o, &mut e);
        (
            c,
            String::from_utf8(o).unwrap(),
            String::from_utf8(e).unwrap(),
        )
    }

    const CAT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../bin/cn_loci.tsv");

    #[test]
    fn marker_regions_lists_every_copy() {
        let (c, o, _) = run(&["marker-regions", "SMN1_27134"]);
        assert_eq!(
            (c, o.as_str()),
            (0, "chr5:70952074-70952074\nchr5:70076654-70076654\n")
        );
        assert_eq!(run(&["marker-regions", "NOPE"]).0, 2);
    }

    #[test]
    fn bed_psv_and_usage() {
        let (c, o, _) = run(&["bed", "--catalog", CAT, "--psv", "--locus", "CYP21A2"]);
        assert_eq!(c, 0);
        assert_eq!(o.lines().count(), 12);
        let (c, _, e) = run(&["bed", "--catalog", CAT, "--psv", "--locus", "HBA"]);
        assert_eq!(
            (c, e.as_str()),
            (2, "ERROR: no psv locus 'HBA' in the catalog\n")
        );
        assert_eq!(run(&["call", "--locus", "HBA"]).0, 2);
    }

    /// A stdout that refuses every write, as a full disk does.
    struct Full;
    impl Write for Full {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("no space left on device"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_region_list_that_cannot_be_written_fails() {
        // A truncated region list would let a step read the marker at one
        // copy, or count over part of a locus; each writer must exit 2.
        for args in [
            &["marker-regions", "SMN1_27134"][..],
            &["bed", "--catalog", CAT, "--locus", "HBA"],
            &["bed", "--catalog", CAT, "--psv", "--locus", "CYP21A2"],
        ] {
            let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
            let mut e = Vec::new();
            let c = main(&args, &mut Full, &mut e);
            let e = String::from_utf8(e).unwrap();
            assert_eq!(c, 2, "{args:?}: {e}");
            assert!(
                e.contains("cannot write the regions: no space left on device"),
                "{e}"
            );
        }
    }
}
