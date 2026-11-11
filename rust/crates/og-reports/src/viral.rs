//! `og viral report` (step 47): the viral DNA in the reads that do not map
//! to the human genome.
//!
//! ```text
//! og viral report --sample S --viral-fasta viral_genomes.fa
//!                 --viral-coverage S_viral_coverage.tsv --viral-depth S_viral_depth.tsv
//!                 --human-coverage S_human_coverage.tsv --out S_viral.tsv
//! ```
//!
//! The unmapped reads are realigned to the viral genomes; `samtools
//! coverage` of that BAM gives reads per virus, `samtools depth -aa` the
//! depth at every viral base, and `samtools coverage` of the sample's BAM
//! the human depth (chr1 to chr22, over the bases with at least one read).
//! One row per viral genome, in the FASTA's order.
//!
//! The call: inherited integrated HHV-6 when the mode of the per-base depth
//! over HHV-6A or HHV-6B is at least 0.30 of the human depth (Liu et al.
//! 2020, PLoS Genet, doi:10.1371/journal.pgen.1008915: one integrated copy
//! per cell gives about 0.5, viraemia 10 to 1000 times less; when both
//! species pass, the higher one is called and the other row reads as
//! reads shared by the two); `lab control (phiX)` for Illumina's spike-in;
//! `viral DNA found` (10 reads or more over 10% of the genome or more);
//! `a few reads`; `not found`. None of these is a diagnosis of infection.
//!
//! The depth file (about 1.6 million lines over the 24 genomes) is read
//! once into a histogram per genome: memory is bounded by the number of
//! distinct depths, not the genome length.
//!
//! What changes on purpose: a number in a samtools table that is not one
//! stops with the file and line (`ERROR:`, exit 1) instead of a traceback,
//! and problems start `ERROR:`, not `viral_report.py:`.

use crate::Error;
use indexmap::IndexMap;
use og_filters::args::Spec;
use og_io::Utf8;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

const MIN_READS: i128 = 10;
const FOUND_BREADTH: f64 = 0.10;
const ICI_RELATIVE_DEPTH: f64 = 0.30;
pub const HEADER: [&str; 12] = [
    "virus",
    "accession",
    "length",
    "reads",
    "breadth",
    "mean_depth",
    "mode_depth",
    "human_depth",
    "relative_depth",
    "copies_per_cell",
    "call",
    "note",
];

/// The call of one viral genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViralCall {
    Ici,
    LabControl,
    Found,
    FewReads,
    NotFound,
}

impl ViralCall {
    fn text(self) -> &'static str {
        match self {
            ViralCall::Ici => "inherited integrated HHV-6 likely",
            ViralCall::LabControl => "lab control (phiX)",
            ViralCall::Found => "viral DNA found",
            ViralCall::FewReads => "a few reads",
            ViralCall::NotFound => "not found",
        }
    }
}

fn hhv6(base: &str) -> Option<&'static str> {
    match base {
        "NC_001664" => Some("HHV-6A"),
        "NC_000898" => Some("HHV-6B"),
        _ => None,
    }
}

const PHIX: &str = "NC_001422";
/// HCV: an RNA virus with no DNA stage.
const RNA_VIRUS: &str = "NC_004102";

fn base(acc: &str) -> &str {
    acc.split('.').next().unwrap_or(acc)
}

/// `[(accession, title)]` of the FASTA, in file order.
fn read_fasta_titles(path: &str) -> Result<Vec<(String, String)>, Error> {
    let mut out = Vec::new();
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(file, Utf8::Replace);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        if let Some(rest) = line.strip_prefix('>') {
            let t = crate::strip(rest);
            let (acc, title) = t.split_once(' ').unwrap_or((t, ""));
            let title = crate::strip(title);
            out.push((
                acc.to_owned(),
                if title.is_empty() { acc } else { title }.to_owned(),
            ));
        }
    }
    if out.is_empty() {
        return Err(Error(format!("{path} has no sequence")));
    }
    Ok(out)
}

/// One `samtools coverage` row.
#[derive(Debug, Clone, Copy)]
struct Cov {
    reads: i128,
    length: i128,
    covbases: i128,
    meandepth: f64,
}

/// `samtools coverage` rows by contig, in file order.
fn read_coverage(path: &str) -> Result<IndexMap<String, Cov>, Error> {
    let mut rows = IndexMap::new();
    let mut header: Option<Vec<String>> = None;
    for (i, line) in crate::read_lines(path, Utf8::Replace)?.iter().enumerate() {
        let cols: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(line)
            .split('\t')
            .collect();
        if line.starts_with('#') {
            header = Some(
                cols.iter()
                    .map(|c| c.trim_start_matches('#').to_owned())
                    .collect(),
            );
            continue;
        }
        let Some(h) = header.as_ref().filter(|h| h.len() == cols.len()) else {
            return Err(Error(format!("{path} is not samtools coverage output")));
        };
        let r: HashMap<&str, &str> = h
            .iter()
            .map(String::as_str)
            .zip(cols.iter().copied())
            .collect();
        let cell = |k: &str| {
            r.get(k).copied().ok_or_else(|| {
                Error(format!(
                    "{path} is not samtools coverage output (no {k} column)"
                ))
            })
        };
        let bad = |k: &str, v: &str| {
            Error(format!(
                "{path} line {}: {k} {} is not a number",
                i + 1,
                og_filters::py_repr(v)
            ))
        };
        let int = |k: &str| -> Result<i128, Error> {
            let v = cell(k)?;
            crate::int(v).ok_or_else(|| bad(k, v))
        };
        let mean = cell("meandepth")?;
        let cov = Cov {
            reads: int("numreads")?,
            length: int("endpos")? - int("startpos")? + 1,
            covbases: int("covbases")?,
            meandepth: crate::float(mean).ok_or_else(|| bad("meandepth", mean))?,
        };
        rows.insert(cell("rname")?.to_owned(), cov);
    }
    Ok(rows)
}

/// chr1 to chr22, with or without the prefix.
fn autosome(name: &str) -> bool {
    let n = name.strip_prefix("chr").unwrap_or(name);
    let b = n.as_bytes();
    match b {
        [d] => (b'1'..=b'9').contains(d),
        [b'1', d] => d.is_ascii_digit(),
        [b'2', d] => (b'0'..=b'2').contains(d),
        _ => false,
    }
}

/// Mean depth over chr1 to chr22's bases with at least one read.
fn human_depth(cov: &IndexMap<String, Cov>) -> Option<f64> {
    let (mut total, mut covered) = (0.0f64, 0.0f64);
    for (name, r) in cov {
        if autosome(name) && r.covbases > 0 {
            total += r.meandepth * r.length as f64;
            covered += r.covbases as f64;
        }
    }
    (covered != 0.0).then(|| total / covered)
}

/// Depths seen over one genome: small depths in a vector, the rest in a map.
#[derive(Default)]
struct Hist {
    small: Vec<u64>,
    big: BTreeMap<i128, u64>,
    positions: u64,
}

const SMALL: usize = 4096;

impl Hist {
    fn add(&mut self, d: i128) {
        self.positions += 1;
        match usize::try_from(d) {
            Ok(i) if i < SMALL => {
                if self.small.len() <= i {
                    self.small.resize(i + 1, 0);
                }
                self.small[i] += 1;
            }
            _ => *self.big.entry(d).or_default() += 1,
        }
    }

    /// (depth, count) of every depth seen, smallest first.
    fn counts(&self) -> Vec<(i128, u64)> {
        let mut v: Vec<(i128, u64)> = self.big.range(..0).map(|(d, c)| (*d, *c)).collect();
        v.extend(
            self.small
                .iter()
                .enumerate()
                .filter(|(_, c)| **c > 0)
                .map(|(d, c)| (d as i128, *c)),
        );
        v.extend(self.big.range(0..).map(|(d, c)| (*d, *c)));
        v
    }
}

/// (breadth, mean, mode) per genome of `lengths`, from `samtools depth -aa`.
fn read_depth(
    path: &str,
    lengths: &IndexMap<String, i128>,
) -> Result<HashMap<String, (f64, f64, i128)>, Error> {
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(file, Utf8::Replace);
    let mut hist: HashMap<String, Hist> = HashMap::new();
    let mut line = String::new();
    let mut n = 0u64;
    let mut last: Option<String> = None;
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        n += 1;
        let mut cols = line.split('\t');
        let (Some(name), Some(_), Some(d), None) =
            (cols.next(), cols.next(), cols.next(), cols.next())
        else {
            return Err(Error(format!(
                "{path} line {n} is not samtools depth output (contig, position, depth)"
            )));
        };
        let d = crate::int(d).ok_or_else(|| {
            Error(format!(
                "{path} line {n}: depth {} is not a whole number",
                og_filters::py_repr(d.trim_end_matches('\n'))
            ))
        })?;
        // The file runs contig by contig: look the contig up only when it changes.
        if last.as_deref() != Some(name) {
            last = Some(name.to_owned());
            hist.entry(name.to_owned()).or_default();
        }
        if let Some(h) = hist.get_mut(name) {
            h.add(d);
        }
    }
    let mut out = HashMap::new();
    for (name, &length) in lengths {
        let seen = hist.get(name).map_or(0, |h| h.positions);
        if seen != 0 && i128::from(seen) != length {
            return Err(Error(format!(
                "{path} has {seen} positions for {name}, its length is {length} (run samtools depth -aa)"
            )));
        }
        if length == 0 {
            return Err(Error(format!("{name} has length 0 in the viral coverage")));
        }
        let counts = match hist.get(name) {
            Some(h) if h.positions > 0 => h.counts(),
            _ => vec![(0, u64::try_from(length).unwrap_or(0))],
        };
        let zero = counts.iter().find(|(d, _)| *d == 0).map_or(0, |(_, c)| *c);
        let breadth = 1.0 - zero as f64 / length as f64;
        let sum: i128 = counts.iter().map(|(d, c)| d * i128::from(*c)).sum();
        let mean = sum as f64 / length as f64;
        // The most common depth; a tie goes to the lower one.
        let mode = counts
            .iter()
            .fold(None::<(i128, u64)>, |best, &(d, c)| match best {
                Some((_, bc)) if bc >= c => best,
                _ => Some((d, c)),
            })
            .map_or(0, |(d, _)| d);
        out.insert(name.clone(), (breadth, mean, mode));
    }
    Ok(out)
}

/// A word character for Python's `\b` on a str.
fn word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// The title without its `, complete genome` tail:
/// `re.sub(r",? (complete (virion )?genome|genome|DNA)\b.*$", "", title)`.
fn short_title(title: &str) -> &str {
    const ALTS: [&str; 4] = ["complete virion genome", "complete genome", "genome", "DNA"];
    let at = |rest: &str| {
        ALTS.iter().any(|a| {
            rest.strip_prefix(a)
                .is_some_and(|after| !after.chars().next().is_some_and(word))
        })
    };
    for (i, _) in title.char_indices() {
        let rest = &title[i..];
        let tail = rest.strip_prefix(", ").or_else(|| rest.strip_prefix(' '));
        if tail.is_some_and(at) {
            return &title[..i];
        }
    }
    title
}

fn short_name(acc: &str, title: &str) -> String {
    if let Some(n) = hhv6(base(acc)) {
        return n.into();
    }
    let t = crate::strip(short_title(title));
    if t.is_empty() { acc } else { t }.to_owned()
}

/// The call and its note for one viral genome (before the HHV-6 pick).
fn call_row(acc: &str, reads: i128, breadth: f64, relative: Option<f64>) -> (ViralCall, String) {
    let b = base(acc);
    if b == PHIX {
        if reads != 0 {
            return (
                ViralCall::LabControl,
                "Illumina adds phiX to the run as a control: these reads were added in the lab."
                    .into(),
            );
        }
        return (
            ViralCall::NotFound,
            "Illumina's phiX control: none of its reads reached this sample's data.".into(),
        );
    }
    if let (Some(_), Some(rel)) = (hhv6(b), relative)
        && rel >= ICI_RELATIVE_DEPTH
    {
        return (
            ViralCall::Ici,
            format!(
                "Depth over the genome at {} of the human depth (Liu et al. 2020 call integration at {} or more; one copy per cell gives about 0.5). Inherited integrated HHV-6 is in every cell and is passed on to about half of the children; it can make a blood test read as an active HHV-6 infection.",
                crate::fixed(rel, 2),
                crate::fixed(ICI_RELATIVE_DEPTH, 2)
            ),
        );
    }
    if reads < MIN_READS {
        return (ViralCall::NotFound, String::new());
    }
    let note = if hhv6(b).is_some() {
        if relative.is_some() {
            "Below the integration threshold of Liu et al. 2020: not the one-copy-per-cell pattern of inherited HHV-6."
        } else {
            "No human depth (no read on chr1 to chr22): integration is not assessed."
        }
    } else if b == RNA_VIRUS {
        "An RNA virus with no DNA stage: DNA reads point to a contamination or an integration (Moustafa et al. 2017)."
    } else {
        ""
    };
    let call = if breadth >= FOUND_BREADTH {
        ViralCall::Found
    } else {
        ViralCall::FewReads
    };
    (call, note.into())
}

/// One row: the cells of [`HEADER`] and the call.
struct Row {
    virus: String,
    cells: Vec<String>,
    call: ViralCall,
    note: String,
    relative: f64,
}

/// The files of `og viral report`.
pub struct Files<'a> {
    pub viral_fasta: &'a str,
    pub viral_coverage: &'a str,
    pub viral_depth: &'a str,
    pub human_coverage: &'a str,
    pub out: &'a str,
}

/// `og viral report`: writes the table; the stdout line.
pub fn report(sample: &str, f: &Files<'_>) -> Result<String, Error> {
    let titles = read_fasta_titles(f.viral_fasta)?;
    let vcov = read_coverage(f.viral_coverage)?;
    let lengths: IndexMap<String, i128> = titles
        .iter()
        .filter_map(|(a, _)| vcov.get(a).map(|c| (a.clone(), c.length)))
        .collect();
    let vdepth = read_depth(f.viral_depth, &lengths)?;
    let hdepth = human_depth(&read_coverage(f.human_coverage)?).filter(|h| *h != 0.0);
    let mut rows = Vec::new();
    for (acc, title) in &titles {
        let Some(c) = vcov.get(acc) else {
            return Err(Error(format!(
                "the viral coverage has no row for {acc}: was the BAM aligned to this FASTA?"
            )));
        };
        let (breadth, mean, mode) = vdepth[acc];
        let relative = hdepth.map(|h| mode as f64 / h);
        let copies = hdepth.map(|h| 2.0 * mean / h);
        let (call, note) = call_row(acc, c.reads, breadth, relative);
        let virus = short_name(acc, title);
        rows.push(Row {
            cells: vec![
                virus.clone(),
                acc.clone(),
                c.length.to_string(),
                c.reads.to_string(),
                crate::fixed(breadth, 4),
                crate::fixed(mean, 3),
                mode.to_string(),
                hdepth.map_or(".".into(), |h| crate::fixed(h, 2)),
                relative.map_or(".".into(), |r| crate::fixed(r, 3)),
                copies.map_or(".".into(), |x| crate::fixed(x, 4)),
            ],
            virus,
            call,
            note,
            relative: relative.unwrap_or(0.0),
        });
    }
    // Both HHV-6 species pass: the higher relative depth is the integrated
    // genome; the other's depth is reads shared by the two.
    let ici: Vec<usize> = (0..rows.len())
        .filter(|&i| rows[i].call == ViralCall::Ici)
        .collect();
    if ici.len() > 1 {
        let keep = ici.iter().copied().fold(ici[0], |k, i| {
            if rows[i].relative > rows[k].relative {
                i
            } else {
                k
            }
        });
        let name = rows[keep].virus.clone();
        for &i in &ici {
            if i != keep {
                rows[i].call = ViralCall::Found;
                rows[i].note = format!(
                    "Reads shared with {name}, which has the higher depth and is the integrated species called."
                );
            }
        }
    }
    let mut w = pycompat::csv::Writer::new(Vec::new(), '\t');
    let io = |e: std::io::Error| Error::io(&e, f.out);
    w.write_row(HEADER).map_err(io)?;
    for r in &rows {
        let mut cells = r.cells.clone();
        cells.push(r.call.text().into());
        cells.push(r.note.clone());
        w.write_row(cells).map_err(io)?;
    }
    crate::write_atomic(f.out, &w.into_inner())?;
    let found: Vec<String> = rows
        .iter()
        .filter(|r| r.call != ViralCall::NotFound)
        .map(|r| format!("{} ({})", r.virus, r.call.text()))
        .collect();
    Ok(format!(
        "{sample}: {} viral genomes; {}\n",
        rows.len(),
        if found.is_empty() {
            "no viral DNA found".to_owned()
        } else {
            found.join("; ")
        }
    ))
}

const SPEC: Spec = Spec {
    prog: "og viral report",
    usage: "usage: og viral report [-h] --sample SAMPLE --viral-fasta VIRAL_FASTA --viral-coverage VIRAL_COVERAGE --viral-depth VIRAL_DEPTH --human-coverage HUMAN_COVERAGE --out OUT\n",
    options: &[
        "--sample",
        "--viral-fasta",
        "--viral-coverage",
        "--viral-depth",
        "--human-coverage",
        "--out",
    ],
    positionals: &[],
    required: &[
        "--sample",
        "--viral-fasta",
        "--viral-coverage",
        "--viral-depth",
        "--human-coverage",
        "--out",
    ],
};

/// `og viral VERB ARGS`: the exit code.
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if verb != Some("report") {
        return crate::unknown_verb("viral", "report", verb, err);
    }
    crate::run(&SPEC, rest, out, err, |a| {
        let files = Files {
            viral_fasta: crate::opt(a, "--viral-fasta"),
            viral_coverage: crate::opt(a, "--viral-coverage"),
            viral_depth: crate::opt(a, "--viral-depth"),
            human_coverage: crate::opt(a, "--human-coverage"),
            out: crate::opt(a, "--out"),
        };
        Ok(report(crate::opt(a, "--sample"), &files).map_err(|e| (e, 1)))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titles_and_names() {
        assert_eq!(
            short_title("Human gammaherpesvirus 4, complete genome"),
            "Human gammaherpesvirus 4"
        );
        assert_eq!(
            short_title("Hepatitis B virus genome sequence"),
            "Hepatitis B virus"
        );
        assert_eq!(
            short_title("Torque teno virus DNA, isolate x"),
            "Torque teno virus"
        );
        assert_eq!(short_title("Some genomes"), "Some genomes");
        assert_eq!(short_name("NC_001664.4", "x"), "HHV-6A");
        assert_eq!(short_name("NC_1", ", complete genome"), "NC_1");
        assert!(autosome("chr22") && autosome("7") && !autosome("chr23") && !autosome("chr0"));
    }

    #[test]
    fn mode_tie_goes_to_the_lower_depth() {
        let d = crate::testutil::Dir::new("viral");
        let p = d.put("d.tsv", "V\t1\t3\nV\t2\t3\nV\t3\t1\nV\t4\t1\n");
        let lengths: IndexMap<String, i128> = [("V".to_owned(), 4)].into_iter().collect();
        let got = read_depth(&p, &lengths).unwrap();
        assert_eq!(got["V"], (1.0, 2.0, 1));
    }

    /// Red first: a relative depth of exactly 0.30 is integrated HHV-6.
    #[test]
    fn red_first_the_ici_threshold_is_inclusive() {
        assert_eq!(
            call_row("NC_000898.1", 0, 0.0, Some(0.30)).0,
            ViralCall::Ici
        );
        assert_eq!(
            call_row("NC_000898.1", 50, 0.5, Some(0.2999)).0,
            ViralCall::Found
        );
        assert_eq!(
            call_row("NC_001422.1", 3, 0.0, None).0,
            ViralCall::LabControl
        );
    }
}
