//! One PGS Catalog scoring file, written as pgsc_calc's custom GRCh38 file.
//!
//! The file is read twice and never held. The first pass reads the header
//! (`#key=value` lines, wherever they are; the last of a key wins), the
//! column line, whether any row is dominant or recessive, and, with a
//! reference, the strand check. The second pass writes the rows and hands
//! their alleles to [`Alleles`]. So a refusal comes before anything is
//! written, as it did when the whole file was read first.

use super::alleles::{Alleles, ScoreAlleles, chrom_rank};
use flate2::Compression;
use flate2::write::GzEncoder;
use og_intervals::fasta::Indexed;
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::text::{isdigit, strip, upper};
use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Columns a score must have.
const REQUIRED: [&str; 4] = ["hm_chr", "hm_pos", "effect_allele", "effect_weight"];
/// Columns that make a score non-additive.
const NON_ADDITIVE: [&str; 3] = ["dosage_0_weight", "dosage_1_weight", "dosage_2_weight"];
/// Columns whose `True` makes a row non-additive, checked in this order.
const MODELS: [&str; 2] = ["is_dominant", "is_recessive"];

/// What happens to a score's strand-ambiguous rows (A/T and C/G pairs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ambiguous {
    /// Written; pgsc_calc drops them (no `--keep-ambiguous`).
    Dropped,
    /// Written and kept: the score passed the strand check.
    Kept,
    /// Left out: the score failed the strand check.
    StrandCheckFailed,
}

impl Ambiguous {
    pub fn as_str(self) -> &'static str {
        match self {
            Ambiguous::Dropped => "dropped",
            Ambiguous::Kept => "kept",
            Ambiguous::StrandCheckFailed => "dropped, strand check failed",
        }
    }
}

/// What formatting one score did.
#[derive(Debug)]
pub struct Summary {
    /// The header's `#pgs_id`, else the file name's id.
    pub pid: String,
    pub written: u64,
    /// Rows dropped for their chromosome: chrY, MT, other contigs, and chrX without `keep_x`.
    pub off_autosomes: u64,
    /// Rows on chrX, written or not.
    pub chrx: u64,
    pub ambiguous: u64,
    pub fate: Ambiguous,
    /// Strand check counts (with a reference only).
    pub agree: u64,
    pub flipped: u64,
}

/// What a score needs from the run.
pub struct Options<'a> {
    pub keep_x: bool,
    /// The label of each id, from the score list.
    pub label: &'a dyn Fn(&str) -> Option<String>,
}

/// `PGS000018` for `PGS000018.txt.gz`, `PGS000018_hmPOS_GRCh38.txt.gz` and the like.
pub fn pgs_id_of(path: &Path) -> String {
    let mut name = path
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    for suffix in [".gz", ".txt", ".tsv"] {
        if let Some(s) = name.strip_suffix(suffix) {
            name = s.to_owned();
        }
    }
    match name.split_once("_hmPOS_") {
        Some((id, _)) => id.to_owned(),
        None => name,
    }
}

/// `str.upper()`, without an allocation for the ASCII it nearly always is.
fn py_upper(s: &str) -> Cow<'_, str> {
    if s.bytes().any(|b| b.is_ascii_lowercase()) || !s.is_ascii() {
        Cow::Owned(if s.is_ascii() {
            s.to_ascii_uppercase()
        } else {
            upper(s)
        })
    } else {
        Cow::Borrowed(s)
    }
}

/// The complement of a single base, as the strand check reads it.
fn complement(b: &str) -> Option<&'static str> {
    match b {
        "A" => Some("T"),
        "T" => Some("A"),
        "C" => Some("G"),
        "G" => Some("C"),
        _ => None,
    }
}

/// A strand-ambiguous SNV: its two alleles are complements (A/T, C/G).
fn ambiguous_pair(ea: &str, oa: &str) -> bool {
    complement(ea).is_some_and(|c| c == oa)
}

/// An other allele the file names as one allele: stripped and upper case;
/// empty for none, `.` or a `/` list.
fn single_allele(a: &str) -> Cow<'_, str> {
    let a = py_upper(strip(a, None));
    if a.is_empty() || a.contains('/') || a == "." {
        Cow::Borrowed("")
    } else {
        a
    }
}

/// Where the columns a row is read by are.
#[derive(Debug, Default)]
struct Columns {
    names: Vec<String>,
    chr: Option<usize>,
    pos: Option<usize>,
    effect: Option<usize>,
    weight: Option<usize>,
    other: Option<usize>,
    inferred: Option<usize>,
    models: [Option<usize>; 2],
}

impl Columns {
    fn new(line: &str) -> Columns {
        let names: Vec<String> = line.split('\t').map(str::to_owned).collect();
        // A name repeated keeps its last column, as a dict built from the line does.
        let at = |k: &str| names.iter().rposition(|n| n == k);
        Columns {
            chr: at("hm_chr"),
            pos: at("hm_pos"),
            effect: at("effect_allele"),
            weight: at("effect_weight"),
            other: at("other_allele"),
            inferred: at("hm_inferOtherAllele"),
            models: [at(MODELS[0]), at(MODELS[1])],
            names,
        }
    }

    fn has(&self, k: &str) -> bool {
        self.names.iter().any(|n| n == k)
    }
}

/// The cells of a row, split once.
struct Row<'a> {
    cells: Vec<&'a str>,
}

impl<'a> Row<'a> {
    fn new(line: &'a str, buf: Vec<&'a str>) -> Row<'a> {
        let mut cells = buf;
        cells.clear();
        cells.extend(line.split('\t'));
        Row { cells }
    }

    /// The stripped cell of a column; empty when the column or the cell is missing.
    fn get(&self, col: Option<usize>) -> &'a str {
        col.and_then(|i| self.cells.get(i))
            .map_or("", |c| strip(c, None))
    }
}

/// A row pgsc_calc gets.
struct Placed<'a> {
    chrom: &'a str,
    rank: u8,
    pos: &'a str,
    effect: Cow<'a, str>,
    other: Cow<'a, str>,
    weight: &'a str,
}

/// Row counts by chromosome, as the row filter makes them.
#[derive(Default)]
struct Counts {
    off: u64,
    x: u64,
}

/// The row as pgsc_calc gets it, or `None` for a row the catalog could not
/// place on GRCh38 or one on a chromosome left out.
fn placed<'a>(r: &Row<'a>, c: &Columns, keep_x: bool, counts: &mut Counts) -> Option<Placed<'a>> {
    let (chrom, pos, weight) = (r.get(c.chr), r.get(c.pos), r.get(c.weight));
    let effect = py_upper(r.get(c.effect));
    if chrom.is_empty() || pos.is_empty() || effect.is_empty() || weight.is_empty() {
        return None;
    }
    let chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    if chrom == "X" {
        counts.x += 1;
    }
    let rank = chrom_rank(chrom).filter(|_| chrom != "X" || keep_x);
    let Some(rank) = rank else {
        counts.off += 1;
        return None;
    };
    let mut other = single_allele(r.get(c.other));
    if other.is_empty() {
        other = single_allele(r.get(c.inferred));
    }
    Some(Placed {
        chrom,
        rank,
        pos,
        effect,
        other,
        weight,
    })
}

/// The strand check: rows whose alleles agree with the reference, and rows
/// that agree only on the other strand.
#[derive(Default)]
struct Strand {
    agree: u64,
    flipped: u64,
    /// Each chromosome name of the score to its name in the FASTA.
    names: HashMap<String, Option<String>>,
}

impl Strand {
    /// Counts one placed row: one whose alleles are two single bases and not
    /// an ambiguous pair, at a position of a sequence the reference has.
    fn count(&mut self, p: &Placed, fasta: &mut Indexed) -> Result<(), String> {
        let (ea, oa) = (p.effect.as_ref(), p.other.as_ref());
        if ea.chars().count() != 1
            || oa.chars().count() != 1
            || ambiguous_pair(ea, oa)
            || !isdigit(p.pos)
        {
            return Ok(());
        }
        let Ok(pos) = u64::try_from(pycompat::num::parse_int(p.pos).unwrap_or(0)) else {
            return Ok(());
        };
        let name = self
            .names
            .entry(p.chrom.to_owned())
            .or_insert_with(|| {
                let alt = format!("chr{}", p.chrom);
                [p.chrom, alt.as_str()]
                    .into_iter()
                    .find(|n| fasta.entry(n).is_some())
                    .map(str::to_owned)
            })
            .as_deref();
        let Some(name) = name else { return Ok(()) };
        if pos == 0 {
            return Ok(());
        }
        let base = fasta
            .base(name, pos - 1)
            .map_err(|e| format!("the reference: {e}"))?;
        let Some(base) = base.map(|b| b.to_ascii_uppercase()) else {
            return Ok(());
        };
        let b = [base];
        let b = std::str::from_utf8(&b).unwrap_or("");
        if complement(b).is_none() {
            return Ok(());
        }
        if b == ea || b == oa {
            self.agree += 1;
        } else if complement(ea) == Some(b) || complement(oa) == Some(b) {
            self.flipped += 1;
        }
        Ok(())
    }
}

/// The lines of a score file, in Python's text mode with `errors="replace"`.
fn lines(path: &Path) -> Result<Lines<Box<dyn std::io::BufRead>>, String> {
    let input = og_io::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Lines::new(input, Utf8::Replace))
}

/// What the first pass learns.
struct Scan {
    header: HashMap<String, String>,
    cols: Option<Columns>,
    models_true: [bool; 2],
    strand: Strand,
    /// A reference read error, reported after the header checks.
    strand_error: Option<String>,
}

/// The first pass.
fn scan(path: &Path, keep_x: bool, mut fasta: Option<&mut Indexed>) -> Result<Scan, String> {
    let mut s = Scan {
        header: HashMap::new(),
        cols: None,
        models_true: [false; 2],
        strand: Strand::default(),
        strand_error: None,
    };
    let mut lines = lines(path)?;
    let mut line = String::new();
    let mut buf: Vec<&'static str> = Vec::new();
    let mut counts = Counts::default();
    while lines
        .read_line(&mut line)
        .map_err(|e| format!("{}: {e}", path.display()))?
    {
        let text = line.strip_suffix('\n').unwrap_or(&line);
        if let Some(h) = text.strip_prefix('#') {
            let (k, v) = h.split_once('=').unwrap_or((h, ""));
            s.header
                .insert(strip(k, None).to_owned(), strip(v, None).to_owned());
            continue;
        }
        let Some(c) = &s.cols else {
            s.cols = Some(Columns::new(text));
            continue;
        };
        if text.is_empty() {
            continue;
        }
        let want_models = c
            .models
            .iter()
            .zip(&s.models_true)
            .any(|(m, t)| m.is_some() && !t);
        if !want_models && fasta.is_none() {
            continue;
        }
        let row = Row::new(text, std::mem::take(&mut buf));
        for (i, m) in c.models.iter().enumerate() {
            if m.is_some() && row.get(*m) == "True" {
                s.models_true[i] = true;
            }
        }
        if let Some(f) = fasta.as_deref_mut()
            && s.strand_error.is_none()
            && let Some(p) = placed(&row, c, keep_x, &mut counts)
        {
            // Reported after the header checks, which come first.
            s.strand_error = s.strand.count(&p, f).err();
        }
        buf = recycle(row.cells);
    }
    Ok(s)
}

/// An emptied cell vector whose allocation the next row reuses.
fn recycle(mut v: Vec<&str>) -> Vec<&'static str> {
    v.clear();
    v.into_iter().map(|_| "").collect()
}

/// Writes `path` as `out_dir/<pid>.txt.gz` and adds its alleles to
/// `alleles`. The error is the message after `ERROR: `.
pub fn format(
    path: &Path,
    out_dir: &Path,
    opts: &Options,
    fasta: Option<&mut Indexed>,
    alleles: &mut Alleles,
) -> Result<Summary, String> {
    let name = path
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    let with_ref = fasta.is_some();
    let mut s = scan(path, opts.keep_x, fasta)?;
    let build = s.header.get("HmPOS_build").map_or("", String::as_str);
    if build != "GRCh38" {
        return Err(format!(
            "{name} has #HmPOS_build='{build}', expected GRCh38 (use the PGS Catalog's Harmonized/<id>_hmPOS_GRCh38 file)"
        ));
    }
    let cols = s.cols.take().unwrap_or_default();
    let missing: Vec<&str> = REQUIRED.into_iter().filter(|k| !cols.has(k)).collect();
    if !missing.is_empty() {
        return Err(format!("{name} has no {} column", missing.join(", ")));
    }
    if NON_ADDITIVE.iter().any(|k| cols.has(k)) {
        return Err(format!(
            "{name} has dosage weights (a non-additive score); only additive scores are supported"
        ));
    }
    if let Some(i) = s.models_true.iter().position(|&t| t) {
        return Err(format!(
            "{name} has {}=True rows (a non-additive score); only additive scores are supported",
            MODELS[i]
        ));
    }
    if let Some(e) = s.strand_error.take() {
        return Err(e);
    }
    let pid = s
        .header
        .get("pgs_id")
        .filter(|p| !p.is_empty())
        .cloned()
        .unwrap_or_else(|| pgs_id_of(path));
    let label = (opts.label)(&pid)
        .or_else(|| {
            s.header
                .get("trait_reported")
                .filter(|t| !t.is_empty())
                .cloned()
        })
        .unwrap_or_else(|| pid.clone());
    let label = label.replace('=', "-").replace('\n', " ");
    let fate = if !with_ref {
        Ambiguous::Dropped
    } else if s.strand.agree > 0 && s.strand.flipped == 0 {
        Ambiguous::Kept
    } else {
        Ambiguous::StrandCheckFailed
    };

    let target = out_dir.join(format!("{pid}.txt.gz"));
    let io_err = |e: std::io::Error| format!("{}: {e}", target.display());
    let file = AtomicFile::create(&target).map_err(io_err)?;
    let mut w = BufWriter::with_capacity(1 << 18, GzEncoder::new(file, Compression::new(1)));
    write!(
        w,
        "#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={label}\n#genome_build=GRCh38\n\
         chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
    )
    .map_err(io_err)?;

    let mut sum = Summary {
        pid,
        written: 0,
        off_autosomes: 0,
        chrx: 0,
        ambiguous: 0,
        fate,
        agree: s.strand.agree,
        flipped: s.strand.flipped,
    };
    let mut counts = Counts::default();
    let mut mine = ScoreAlleles::default();
    let mut lines = lines(path)?;
    let mut line = String::new();
    let mut seen_cols = false;
    let mut buf: Vec<&'static str> = Vec::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| format!("{}: {e}", path.display()))?
    {
        let text = line.strip_suffix('\n').unwrap_or(&line);
        if text.starts_with('#') {
            continue;
        }
        if !seen_cols {
            seen_cols = true;
            continue;
        }
        if text.is_empty() {
            continue;
        }
        let row = Row::new(text, std::mem::take(&mut buf));
        if let Some(p) = placed(&row, &cols, opts.keep_x, &mut counts) {
            let amb = ambiguous_pair(&p.effect, &p.other);
            if amb {
                sum.ambiguous += 1;
            }
            if !(amb && fate == Ambiguous::StrandCheckFailed) {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}",
                    p.chrom, p.pos, p.effect, p.other, p.weight
                )
                .map_err(io_err)?;
                sum.written += 1;
                alleles.push(&mut mine, p.rank, p.pos, &p.effect);
                if !p.other.is_empty() {
                    alleles.push(&mut mine, p.rank, p.pos, &p.other);
                }
            }
        }
        buf = recycle(row.cells);
    }
    sum.off_autosomes = counts.off;
    sum.chrx = counts.x;
    if sum.written == 0 {
        // The temporary file goes with the dropped writer.
        return Err(format!(
            "{name} has no autosomal row with hm_chr, hm_pos, effect_allele and effect_weight"
        ));
    }
    let gz = w.into_inner().map_err(|e| io_err(e.into_error()))?;
    gz.finish().and_then(AtomicFile::commit).map_err(io_err)?;
    alleles.add(mine);
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_from_file_names() {
        for (f, id) in [
            ("PGS000018.txt.gz", "PGS000018"),
            ("PGS000018_hmPOS_GRCh38.txt.gz", "PGS000018"),
            ("PGS1.tsv", "PGS1"),
            ("x/PGS2.txt", "PGS2"),
        ] {
            assert_eq!(pgs_id_of(Path::new(f)), id);
        }
    }

    #[test]
    fn single_alleles() {
        assert_eq!(single_allele(" g "), "G");
        assert_eq!(single_allele("A/G"), "");
        assert_eq!(single_allele("."), "");
        assert_eq!(single_allele("ß"), "SS");
    }

    #[test]
    fn ambiguous_pairs() {
        assert!(ambiguous_pair("A", "T"));
        assert!(ambiguous_pair("G", "C"));
        assert!(!ambiguous_pair("A", "G"));
        assert!(!ambiguous_pair("AT", "TA"));
    }
}
