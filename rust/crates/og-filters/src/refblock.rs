//! `og pgx refblock-fill [--min-gq 20] POSITIONS GVCF < expanded.vcf`:
//! PharmCAT's multi-base positions called 0/0 from the gVCF's reference
//! blocks (replaces `bin/pgx_refblock_fill.py`).
//!
//! Step 07 expands the gVCF's reference blocks over PharmCAT's regions
//! with `bcftools convert --gvcf2vcf`, one record per base. At a PharmCAT
//! position whose REF or ALT is longer than one base PharmCAT sets such a
//! record aside, so those alleles are uncallable even where the gVCF shows
//! confident reference over the whole span. This filter writes PharmCAT's
//! own record as 0/0 (INFO flag `PGX_REFBLOCK`) in place of the one-base
//! record at that position, when every base of the REF span lies in
//! confident reference blocks (ALT only `<*>`, `<X>`, `<NON_REF>` or `.`;
//! GT all 0; FILTER PASS or `.`; GQ `--min-gq` or more) and no other gVCF
//! record touches the span. GT is the block's at the first base, GQ the
//! lowest over the span. A summary goes to stderr; exit 2 on unreadable
//! input.
//!
//! Fixes against the Python, each with a golden in `tests/parity`:
//! - when the expanded stream had a record at a filled position that was
//!   not a one-base reference block, the Python wrote that record and then
//!   its 0/0 record as well: two records at one position. The gVCF's record
//!   is kept, the fill dropped, and a line says so;
//! - two PharmCAT positions with the same start shared one coverage count,
//!   so a short one could fail on the longer one's bases, or be filled on
//!   the longer one's success. Each position now has its own;
//! - a record after a last line without a newline was glued onto it;
//! - a malformed record on stdin, or a gVCF record without its 8 columns,
//!   was a traceback (exit 1); now a message naming the line, exit 2.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::io::{self, BufRead, Write};

use og_io::{Lines, Utf8};

use crate::args::{self, ArgError, Spec};

const SPEC: Spec = Spec {
    prog: "og pgx refblock-fill",
    usage: "usage: og pgx refblock-fill [--min-gq N] POSITIONS GVCF < expanded.vcf > out.vcf\n",
    options: &["--min-gq"],
    positionals: &["positions", "gvcf"],
    required: &[],
};

pub const MIN_GQ: i128 = 20;
const BLOCK_ALTS: [&str; 4] = ["<*>", "<X>", "<NON_REF>", "."];
const FLAG: &str = "PGX_REFBLOCK";
/// The header line of the flag. The script name in it is kept, so the file
/// PharmCAT reads does not change.
const FLAG_HEADER: &str = "##INFO=<ID=PGX_REFBLOCK,Number=0,Type=Flag,Description=\"0/0 written by pgx_refblock_fill.py: reference blocks cover the whole REF span\">\n";
const PREFIX: &str = "pgx_refblock_fill";

/// Why the inputs could not be read.
#[derive(Debug)]
pub enum RefblockError {
    Io(String, io::Error),
    /// A line the filter cannot read: file, line (1-based) and why.
    Line(String, u64, String),
    /// The gVCF's sample columns, when not one.
    Samples(usize),
    NoHeader,
}

impl fmt::Display for RefblockError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RefblockError::Io(path, e) => write!(f, "{}", crate::os_error(e, path)),
            RefblockError::Line(path, line, why) => write!(f, "{path}:{line}: {why}"),
            RefblockError::Samples(n) => write!(f, "the gVCF must hold one sample, it holds {n}"),
            RefblockError::NoHeader => f.write_str("the gVCF has no #CHROM header line"),
        }
    }
}

/// One PharmCAT position whose REF or an ALT is longer than one base.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Target {
    pub start: i128,
    /// The last base of REF.
    pub end: i128,
    pub id: String,
    pub ref_: String,
    pub alt: String,
    pub gene: String,
}

/// The positions per chromosome (in the order first seen), each list
/// sorted, and how many there are.
pub struct Targets {
    chroms: Vec<(String, Vec<Target>)>,
    index: HashMap<String, usize>,
    pub count: usize,
}

impl Targets {
    fn of(&self, chrom: &str) -> Option<&[Target]> {
        self.index.get(chrom).map(|&i| self.chroms[i].1.as_slice())
    }
}

fn lines_of(path: &str) -> Result<Lines<Box<dyn BufRead>>, RefblockError> {
    og_io::open(path)
        .map(|r| Lines::new(r, Utf8::Strict))
        .map_err(|e| RefblockError::Io(path.to_owned(), e))
}

fn int_at(path: &str, line: u64, what: &str, s: &str) -> Result<i128, RefblockError> {
    pycompat::num::parse_int(s).map_err(|_| {
        RefblockError::Line(
            path.to_owned(),
            line,
            format!("{what} {} is not a whole number", crate::py_repr(s)),
        )
    })
}

/// Reads PharmCAT's positions VCF (plain or bgzip).
pub fn read_positions(path: &str) -> Result<Targets, RefblockError> {
    let mut t = Targets {
        chroms: Vec::new(),
        index: HashMap::new(),
        count: 0,
    };
    let mut number = 0;
    for line in lines_of(path)? {
        let line = line.map_err(|e| RefblockError::Io(path.to_owned(), e))?;
        number += 1;
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            continue;
        }
        let ref_ = c[3];
        if crate::chars(ref_) == 1 && c[4].split(',').all(|a| crate::chars(a) == 1) {
            continue;
        }
        let gene = c
            .get(7)
            .and_then(|info| info.split(';').find_map(|kv| kv.strip_prefix("PX=")))
            .unwrap_or_default();
        let start = int_at(path, number, "POS", c[1])?;
        let len = i128::try_from(crate::chars(ref_)).unwrap_or(i128::MAX);
        let target = Target {
            start,
            end: start + len - 1,
            id: c[2].to_owned(),
            ref_: ref_.to_owned(),
            alt: c[4].to_owned(),
            gene: gene.to_owned(),
        };
        let i = *t.index.entry(c[0].to_owned()).or_insert_with(|| {
            t.chroms.push((c[0].to_owned(), Vec::new()));
            t.chroms.len() - 1
        });
        t.chroms[i].1.push(target);
        t.count += 1;
    }
    for (_, v) in &mut t.chroms {
        v.sort();
    }
    Ok(t)
}

/// Why a gVCF record does or does not count as confident reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Block {
    /// A confident reference block (or one-base reference record).
    Confident {
        gq: i128,
    },
    NotConfident,
}

/// A confident reference block: ALT only block alleles, GT all 0, FILTER
/// PASS or `.`, GQ `min_gq` or more.
pub fn block(alt: &str, filter: &str, gt: &str, gq: Option<i128>, min_gq: i128) -> Block {
    let alts_ok = alt.split(',').all(|a| BLOCK_ALTS.contains(&a));
    let gt_ok = gt.split(['/', '|']).all(|a| a == "0");
    match gq {
        Some(gq) if alts_ok && gt_ok && (filter == "." || filter == "PASS") && gq >= min_gq => {
            Block::Confident { gq }
        }
        _ => Block::NotConfident,
    }
}

/// What the gVCF shows over one target.
#[derive(Debug, Clone)]
struct Coverage {
    /// One bit per base of the REF span.
    covered: Vec<u64>,
    bad: bool,
    gt_at_start: Option<String>,
    min_gq: Option<i128>,
}

impl Coverage {
    fn new(span: usize) -> Self {
        Coverage {
            covered: vec![0; span.div_ceil(64)],
            bad: false,
            gt_at_start: None,
            min_gq: None,
        }
    }

    fn cover(&mut self, from: usize, to: usize) {
        for k in from..=to {
            self.covered[k / 64] |= 1 << (k % 64);
        }
    }

    fn count(&self) -> usize {
        self.covered.iter().map(|w| w.count_ones() as usize).sum()
    }
}

/// `(?:^|;)END=(\d+)` searched in INFO.
fn info_end(info: &str) -> Option<&str> {
    let b = info.as_bytes();
    let mut from = 0;
    while let Some(k) = info[from..].find("END=") {
        let at = from + k;
        if at == 0 || b[at - 1] == b';' {
            let digits = b[at + 4..]
                .iter()
                .take_while(|c| c.is_ascii_digit())
                .count();
            if digits > 0 {
                return Some(&info[at + 4..at + 4 + digits]);
            }
        }
        from = at + 1;
    }
    None
}

/// The gVCF's contig order and, per target that the gVCF shows reference
/// over, its fill: `(chrom index, target index) -> (GT, lowest GQ)`.
type Filled = HashMap<(usize, usize), (String, i128)>;

fn scan_gvcf(
    path: &str,
    targets: &Targets,
    min_gq: i128,
) -> Result<(Vec<String>, Filled), RefblockError> {
    let mut order: Vec<String> = Vec::new();
    let mut known: HashSet<String> = HashSet::new();
    let mut state: HashMap<(usize, usize), Coverage> = HashMap::new();
    let mut samples = false;
    let mut number = 0;
    let line_err = |line: u64, why: String| RefblockError::Line(path.to_owned(), line, why);
    for line in lines_of(path)? {
        let line = line.map_err(|e| RefblockError::Io(path.to_owned(), e))?;
        number += 1;
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let id: String = rest.chars().take_while(|&c| c != ',' && c != '>').collect();
            if !id.is_empty() && known.insert(id.clone()) {
                order.push(id);
            }
            continue;
        }
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            let n = line.split('\t').count().saturating_sub(9);
            if n != 1 {
                return Err(RefblockError::Samples(n));
            }
            samples = true;
            continue;
        }
        if !samples {
            return Err(RefblockError::NoHeader);
        }
        let c: Vec<&str> = line.split('\t').collect();
        let chrom = c[0];
        if known.insert(chrom.to_owned()) {
            order.push(chrom.to_owned());
        }
        let Some(ci) = targets.index.get(chrom).copied() else {
            continue;
        };
        if c.len() < 8 {
            return Err(line_err(
                number,
                format!("a record with {} columns, not 8 or more", c.len()),
            ));
        }
        let pos = int_at(path, number, "POS", c[1])?;
        let end = match info_end(c[7]) {
            Some(e) => int_at(path, number, "END", e)?,
            None => pos + i128::try_from(crate::chars(c[3])).unwrap_or(i128::MAX) - 1,
        };
        let ts = targets.of(chrom).unwrap_or_default();
        // Targets are sorted by start: those starting at or before END,
        // and ending at or after POS, overlap the record.
        let upto = ts.partition_point(|t| t.start <= end);
        let hits: Vec<usize> = (0..upto).filter(|&j| ts[j].end >= pos).collect();
        if hits.is_empty() {
            continue;
        }
        let mut fmt: HashMap<&str, &str> = HashMap::new();
        if let (Some(keys), Some(vals)) = (c.get(8), c.get(9)) {
            for (k, v) in keys.split(':').zip(vals.split(':')) {
                fmt.insert(k, v);
            }
        }
        let gt = fmt.get("GT").copied().unwrap_or(".");
        let gq = fmt.get("GQ").and_then(|g| pycompat::num::parse_int(g).ok());
        let kind = block(c[4], c[6], gt, gq, min_gq);
        for j in hits {
            let t = &ts[j];
            let span = usize::try_from(t.end - t.start + 1).unwrap_or(0);
            let s = state.entry((ci, j)).or_insert_with(|| Coverage::new(span));
            let Block::Confident { gq } = kind else {
                s.bad = true;
                continue;
            };
            let (lo, hi) = (pos.max(t.start), end.min(t.end));
            if lo <= hi {
                // Both inside the span: lo >= start and hi <= end.
                let off = |x: i128| usize::try_from(x - t.start).unwrap_or(0);
                s.cover(off(lo), off(hi));
            }
            if lo == t.start {
                s.gt_at_start = Some(gt.to_owned());
            }
            s.min_gq = Some(s.min_gq.map_or(gq, |m| m.min(gq)));
        }
    }
    if !samples {
        return Err(RefblockError::NoHeader);
    }
    let mut ok = Filled::new();
    for ((ci, j), s) in state {
        let t = &targets.chroms[ci].1[j];
        let span = usize::try_from(t.end - t.start + 1).unwrap_or(0);
        let covered = s.count();
        if let (false, Some(gt), Some(gq)) = (s.bad, s.gt_at_start, s.min_gq)
            && covered == span
            && !gt.is_empty()
        {
            ok.insert((ci, j), (gt, gq));
        }
    }
    Ok((order, ok))
}

/// A 0/0 record to write: its place in the stream, the line and its label
/// in the summary.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Fill {
    key: (usize, i128),
    line: String,
    label: String,
}

/// What [`run`] wrote: the labels of the records it filled, and the
/// positions where it kept the stream's own record instead.
#[derive(Debug, Default)]
pub struct Outcome {
    pub filled: Vec<String>,
    pub kept: Vec<String>,
}

/// Fills the expanded VCF from `input` into `out`.
pub fn run(
    positions: &str,
    gvcf: &str,
    min_gq: i128,
    input: impl BufRead,
    mut out: impl Write,
) -> Result<(usize, Outcome), RefblockError> {
    let targets = read_positions(positions)?;
    let (order, ok) = scan_gvcf(gvcf, &targets, min_gq)?;
    let rank: HashMap<&str, usize> = order
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();
    let mut fills = Vec::new();
    for (ci, (chrom, ts)) in targets.chroms.iter().enumerate() {
        for (j, t) in ts.iter().enumerate() {
            if let Some((gt, gq)) = ok.get(&(ci, j)) {
                fills.push(Fill {
                    key: (
                        rank.get(chrom.as_str()).copied().unwrap_or(rank.len()),
                        t.start,
                    ),
                    line: format!(
                        "{chrom}\t{}\t{}\t{}\t{}\t.\tPASS\t{FLAG}\tGT:GQ\t{gt}:{gq}\n",
                        t.start, t.id, t.ref_, t.alt
                    ),
                    label: format!(
                        "{} {chrom}:{}",
                        if t.gene.is_empty() { "?" } else { &t.gene },
                        t.start
                    ),
                });
            }
        }
    }
    fills.sort();

    let stdin_name = "<stdin>";
    let w = |e| RefblockError::Io("<stdout>".to_owned(), e);
    let mut outcome = Outcome::default();
    // The key of the last record written as it was: a fill at that key
    // would be a second record at one position.
    let mut kept_key = None;
    // Whether the last line written ended without a newline.
    let mut open_line = false;
    let emit = |f: &Fill,
                kept_key: Option<(usize, i128)>,
                open_line: &mut bool,
                out: &mut dyn Write,
                o: &mut Outcome|
     -> io::Result<()> {
        if kept_key == Some(f.key) {
            o.kept.push(
                f.label
                    .split_once(' ')
                    .map_or(&*f.label, |(_, p)| p)
                    .to_owned(),
            );
            return Ok(());
        }
        if *open_line {
            out.write_all(b"\n")?;
            *open_line = false;
        }
        out.write_all(f.line.as_bytes())?;
        o.filled.push(f.label.clone());
        Ok(())
    };
    let mut lines = crate::TextLines::stdin(input);
    let mut line = String::new();
    let mut i = 0;
    let mut number = 0;
    while lines
        .read_line(&mut line)
        .map_err(|e| RefblockError::Io(stdin_name.to_owned(), e))?
    {
        number += 1;
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                out.write_all(FLAG_HEADER.as_bytes()).map_err(w)?;
            }
            out.write_all(line.as_bytes()).map_err(w)?;
            open_line = !line.ends_with('\n');
            continue;
        }
        let c: Vec<&str> = line.splitn(6, '\t').collect();
        let short = |need: usize| {
            RefblockError::Line(
                stdin_name.to_owned(),
                number,
                format!("a record with {} columns, not {need} or more", c.len()),
            )
        };
        if c.len() < 2 {
            return Err(short(2));
        }
        let key = (
            rank.get(c[0]).copied().unwrap_or(rank.len()),
            int_at(stdin_name, number, "POS", c[1])?,
        );
        while i < fills.len() && fills[i].key < key {
            emit(&fills[i], kept_key, &mut open_line, &mut out, &mut outcome).map_err(w)?;
            i += 1;
        }
        if i < fills.len() && fills[i].key == key && c.len() < 5 {
            return Err(short(5));
        }
        if i < fills.len()
            && fills[i].key == key
            && crate::chars(c[3]) == 1
            && BLOCK_ALTS.contains(&c[4])
        {
            emit(&fills[i], None, &mut open_line, &mut out, &mut outcome).map_err(w)?;
            i += 1;
            continue;
        }
        out.write_all(line.as_bytes()).map_err(w)?;
        open_line = !line.ends_with('\n');
        kept_key = Some(key);
    }
    for f in &fills[i..] {
        emit(f, kept_key, &mut open_line, &mut out, &mut outcome).map_err(w)?;
    }
    out.flush().map_err(w)?;
    Ok((targets.count, outcome))
}

/// `og pgx refblock-fill`: exit 0, or 2 for a usage error or unreadable
/// input.
pub fn cli(argv: &[String]) -> u8 {
    let parse = || -> Result<(String, String, i128), ArgError> {
        let a = args::parse(&SPEC, argv)?;
        let min_gq = a.typed("--min-gq", MIN_GQ, "int", args::int)?;
        Ok((a.pos[0].clone(), a.pos[1].clone(), min_gq))
    };
    let (positions, gvcf, min_gq) = match parse() {
        Ok(p) => p,
        Err(e) => return args::report(&SPEC, &e),
    };
    let stdin = io::stdin().lock();
    let stdout = io::BufWriter::with_capacity(1 << 16, io::stdout().lock());
    match run(&positions, &gvcf, min_gq, stdin, stdout) {
        Ok((n, o)) => {
            for p in &o.kept {
                eprintln!(
                    "{PREFIX}: kept the gVCF's record at {p} (not a reference block), not a 0/0 record beside it"
                );
            }
            let list = if o.filled.is_empty() {
                String::new()
            } else {
                format!(": {}", o.filled.join(", "))
            };
            eprintln!(
                "{PREFIX}: {} of {n} PharmCAT positions longer than one base called 0/0 from reference blocks with GQ >= {min_gq}{list}; the others stay as the gVCF has them (a variant call, or missing)",
                o.filled.len()
            );
            0
        }
        Err(e) => {
            eprintln!("{PREFIX}: ERROR: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn red_first_a_block_at_the_floor_counts() {
        assert_eq!(
            block("<*>", "PASS", "0/0", Some(20), 20),
            Block::Confident { gq: 20 }
        );
        assert_eq!(
            block("<*>", "PASS", "0/0", Some(19), 20),
            Block::NotConfident
        );
        assert_eq!(
            block("<NON_REF>,.", ".", "0", Some(30), 20),
            Block::Confident { gq: 30 }
        );
        assert_eq!(
            block("<*>", "PASS", "./.", Some(30), 20),
            Block::NotConfident
        );
        assert_eq!(
            block("<*>", "LowGQ", "0/0", Some(30), 20),
            Block::NotConfident
        );
        assert_eq!(
            block("T,<*>", "PASS", "0/0", Some(30), 20),
            Block::NotConfident
        );
        assert_eq!(block("<*>", "PASS", "0|0", None, 20), Block::NotConfident);
    }

    #[test]
    fn end_is_read_as_the_python_regex_reads_it() {
        assert_eq!(info_end("END=12;X=1"), Some("12"));
        assert_eq!(info_end("XEND=5;END=7x"), Some("7"));
        assert_eq!(info_end("END=;END=9"), Some("9"));
        assert_eq!(info_end("BLOCKEND=4"), None);
    }

    #[test]
    fn coverage_counts_each_base_once() {
        let mut c = Coverage::new(70);
        c.cover(0, 9);
        c.cover(5, 69);
        assert_eq!(c.count(), 70);
    }
}
