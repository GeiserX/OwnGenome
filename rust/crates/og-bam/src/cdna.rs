//! Spliced cDNA in the reads: `og cdna check` and `og cdna scan`.
//!
//! Spliced cDNA or expression-vector DNA (a known lab contaminant), or a
//! processed retrocopy, has no introns. A read across one of its exon
//! junctions aligns as two pieces: the first clipped at an exon's end, the
//! rest (a supplementary alignment, named in the SA tag) starting at a later
//! exon's start. A spliced aligner writes the same read as one alignment
//! with an N gap. Genomic reads run on into the intron.
//!
//! `og cdna check` (was `bin/clip_check.py`) asks, at each ClinVar or CPSR
//! P/LP hit, whether the alt reads look like that: it counts the alt and
//! ref reads clipped at least [`MIN_CLIP`] bases within [`EDGE_SLACK`] of an
//! exon edge on the same side, in the alignment or one of its SA pieces,
//! and flags the hit when at least [`MIN_ALT_READS`] reads carry the alt
//! allele, at least [`ALT_SHARE_MIN`] of them are clipped and at most
//! [`REF_SHARE_MAX`] of the ref reads are. With `--construct` a hit at a
//! known construct site is named, whatever the verdict.
//!
//! `og cdna scan` (was `bin/cdna_scan.py`) counts the same signature over
//! every exon of the genome: a read counts for a junction (exon end to a
//! later exon start, [`MIN_INTRON`] to [`MAX_INTRON`] bases on) from two of
//! its pieces on one strand in read order, or from an N gap. Junctions with
//! [`MIN_READS`] reads within [`LOCUS_GAP`] of each other form a locus,
//! listed with [`MIN_JUNCTIONS`] or more (ten Berk de Boer 2023, Front
//! Genet, PMID 37323659). Each record is read once, the contigs on several
//! threads, and the loci are listed in the BAM's contig order, then by
//! position.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Write};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use og_io::{AtomicFile, Lines, Utf8};
use og_tsv::{Options, Table};
use pycompat::num::parse_int;
use pycompat::text::{is_decimal, is_space, isdigit, lower, strip, upper};
use rust_htslib::bam::{self, FetchDefinition, Read as _, Record};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::args::Args;
use crate::reads::{
    Filter, ReadView, Seen, Site, allele_at, complex_allele, flag, op, ops, sa_pieces, segment,
};

/// The shortest clip, in bases, that counts.
pub const MIN_CLIP: i64 = 10;
/// How far a clip point may lie from an exon edge, in bases.
pub const EDGE_SLACK: i64 = 3;
/// The fewest alt reads a hit is judged on.
pub const MIN_ALT_READS: u64 = 2;
/// The share of alt reads clipped at an edge that flags a hit.
pub const ALT_SHARE_MIN: f64 = 0.5;
/// The largest share of ref reads clipped there that still flags it.
pub const REF_SHARE_MAX: f64 = 0.2;
/// The note of a flagged hit.
pub const FLAG_NOTE: &str = "possible cDNA/transgene contamination: check the reads";
/// The columns of the check's table.
pub const CHECK_COLUMNS: [&str; 11] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "alt_reads",
    "alt_clipped",
    "ref_reads",
    "ref_clipped",
    "cdna_check",
    "cdna_note",
    "construct_site",
];

/// Reads a junction needs.
pub const MIN_READS: usize = 2;
/// Junctions a locus needs.
pub const MIN_JUNCTIONS: usize = 2;
/// The shortest intron of a junction.
pub const MIN_INTRON: i64 = 30;
/// The longest intron of a junction.
pub const MAX_INTRON: i64 = 1_000_000;
/// Junctions this close (donor after the last acceptor) form one locus.
pub const LOCUS_GAP: i64 = 20_000;
/// The columns of the scan's table.
pub const SCAN_COLUMNS: [&str; 6] = [
    "chrom",
    "start",
    "end",
    "junctions",
    "reads",
    "junction_list",
];

// --- exon edges ---------------------------------------------------------------------

/// Exon edges per contig (sorted, each once) and the exons merged, from a
/// BED of 0-based, half-open intervals. Lines without two digit fields
/// after the contig are skipped.
#[derive(Debug, Default)]
pub struct Exons {
    starts: HashMap<String, Vec<i64>>,
    ends: HashMap<String, Vec<i64>>,
    regions: HashMap<String, Vec<(i64, i64)>>,
}

fn int64(s: &str) -> Result<i64, String> {
    parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| format!("invalid literal for int() with base 10: {s:?}"))
}

impl Exons {
    pub fn read(path: &str) -> Result<Exons, String> {
        let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
        let mut spans: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for line in Lines::new(BufReader::new(f), Utf8::Strict) {
            let line = line.map_err(|e| format!("{path}: {e}"))?;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 || !isdigit(f[1]) || !isdigit(f[2]) {
                continue;
            }
            spans
                .entry(f[0].to_owned())
                .or_default()
                .push((int64(f[1])?, int64(f[2])?));
        }
        if spans.is_empty() {
            return Err(format!("no exon intervals in {path}"));
        }
        let mut ex = Exons::default();
        for (contig, mut v) in spans {
            let sorted_unique = |mut e: Vec<i64>| {
                e.sort_unstable();
                e.dedup();
                e
            };
            ex.starts.insert(
                contig.clone(),
                sorted_unique(v.iter().map(|s| s.0).collect()),
            );
            ex.ends.insert(
                contig.clone(),
                sorted_unique(v.iter().map(|s| s.1).collect()),
            );
            v.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::new();
            for (a, b) in v {
                match merged.last_mut() {
                    Some(m) if a <= m.1 => m.1 = m.1.max(b),
                    _ => merged.push((a, b)),
                }
            }
            ex.regions.insert(contig, merged);
        }
        Ok(ex)
    }

    fn starts(&self, contig: &str) -> &[i64] {
        self.starts.get(contig).map_or(&[], Vec::as_slice)
    }

    fn ends(&self, contig: &str) -> &[i64] {
        self.ends.get(contig).map_or(&[], Vec::as_slice)
    }
}

/// The edge within [`EDGE_SLACK`] of `x`, the nearest (the lower one on a tie).
pub fn nearest(edges: &[i64], x: i64) -> Option<i64> {
    let i = edges.partition_point(|&e| e < x - EDGE_SLACK);
    edges[i..]
        .iter()
        .take_while(|&&e| e <= x + EDGE_SLACK)
        .fold(None, |best: Option<i64>, &e| match best {
            Some(b) if (b - x).abs() <= (e - x).abs() => Some(b),
            _ => Some(e),
        })
}

// --- the check ----------------------------------------------------------------------

/// Whether an alignment at `start` with `cigar` on `contig` has a clip of
/// [`MIN_CLIP`] or more at an exon edge on the same side: a left clip at an
/// exon start, a right clip at an exon end.
fn edge_clipped(ex: &Exons, contig: &str, start: i64, cigar: &[(u8, i64)]) -> bool {
    let s = segment(start, cigar.iter().copied());
    (s.left_clip >= MIN_CLIP && nearest(ex.starts(contig), s.start).is_some())
        || (s.right_clip >= MIN_CLIP && nearest(ex.ends(contig), s.end).is_some())
}

/// Whether the read or one of its SA pieces is clipped at an exon edge.
fn read_clipped(ex: &Exons, contig: &str, r: &impl ReadView) -> Result<bool, String> {
    let own: Vec<(u8, i64)> = ops(r.raw_cigar()).collect();
    if edge_clipped(ex, contig, r.pos(), &own) {
        return Ok(true);
    }
    for p in sa_pieces(r.sa().unwrap_or_default()) {
        let p = p?;
        if edge_clipped(ex, p.contig, p.pos0, &p.cigar) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The four counts of one hit: alt reads, alt clipped, ref reads, ref clipped.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    pub alt: u64,
    pub alt_clipped: u64,
    pub reference: u64,
    pub ref_clipped: u64,
}

/// Counts the reads at a hit. A read counts once (its name and mate
/// number), as its first alignment that shows either allele; supplementary
/// alignments count, since a cDNA piece is often one.
pub fn count_reads<R: ReadView>(
    reads: impl IntoIterator<Item = Result<R, String>>,
    ex: &Exons,
    contig: &str,
    site: &Site,
) -> Result<Counts, String> {
    let mut seen: FxHashSet<(Vec<u8>, u16)> = FxHashSet::default();
    let mut n = Counts::default();
    for r in reads {
        let r = r?;
        if !Filter::CDNA.keeps(&r) {
            continue;
        }
        let key = (r.qname().to_vec(), r.flags() & (flag::READ1 | flag::READ2));
        if seen.contains(&key) {
            continue;
        }
        let Some(a) = allele_at(&r, site) else {
            continue;
        };
        let clipped = read_clipped(ex, contig, &r)?;
        seen.insert(key);
        let (reads, clip) = match a {
            Seen::Alt => (&mut n.alt, &mut n.alt_clipped),
            Seen::Ref => (&mut n.reference, &mut n.ref_clipped),
        };
        *reads += 1;
        *clip += u64::from(clipped);
    }
    Ok(n)
}

/// (cdna_check, cdna_note) from the counts.
pub fn verdict(n: Counts) -> (&'static str, String) {
    if n.alt < MIN_ALT_READS {
        return (
            "not_checked",
            format!("fewer than {MIN_ALT_READS} reads carry the alt allele"),
        );
    }
    #[allow(clippy::cast_precision_loss)] // read counts are far below 2^52
    let share = |a: u64, b: u64, s: f64| (a as f64, s * b as f64);
    let (ac, an) = share(n.alt_clipped, n.alt, ALT_SHARE_MIN);
    let (rc, rn) = share(n.ref_clipped, n.reference, REF_SHARE_MAX);
    if ac >= an && (n.reference == 0 || rc <= rn) {
        ("flagged", FLAG_NOTE.to_owned())
    } else {
        ("clear", String::new())
    }
}

/// One hit: the contig name as the table gives it, the position as Python's
/// `str(int(pos))`, REF and ALT as given.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hit {
    pub chrom: String,
    pub pos: i64,
    pub reference: String,
    pub alt: String,
}

fn cell<'a>(row: &'a og_tsv::Row, col: &str) -> Result<&'a str, String> {
    row.text(col).map_err(|e| e.to_string())
}

/// Step 06's hits table: the columns chrom, pos, ref and alt, by name.
pub fn read_hits(path: &str) -> Result<Vec<Hit>, String> {
    let mut t = Table::open(path).map_err(|e| e.to_string())?;
    if let Some(h) = t.header().map_err(|e| e.to_string())?
        && let Some(miss) = ["chrom", "pos", "ref", "alt"]
            .iter()
            .find(|c| !h.iter().any(|x| x == *c))
    {
        return Err(format!("{path} has no column {miss}"));
    }
    t.map(|row| {
        let row = row.map_err(|e| e.to_string())?;
        Ok(Hit {
            chrom: cell(&row, "chrom")?.to_owned(),
            pos: int64(cell(&row, "pos")?)?,
            reference: cell(&row, "ref")?.to_owned(),
            alt: cell(&row, "alt")?.to_owned(),
        })
    })
    .collect()
}

/// A CPSR or ClinVar class in one spelling: `Likely_Pathogenic` and `LP`
/// read `likely pathogenic`; an empty or NA cell reads empty.
pub fn norm_class(v: Option<&str>) -> String {
    let v = lower(strip(&v.unwrap_or_default().replace('_', " "), None));
    match v.as_str() {
        "p" => "pathogenic",
        "lp" => "likely pathogenic",
        "b" => "benign",
        "lb" => "likely benign",
        "na" | "." => "",
        other => other,
    }
    .to_owned()
}

const PLP: [&str; 2] = ["pathogenic", "likely pathogenic"];

/// `chrom:g.posREF>ALT`, as CPSR's GENOMIC_CHANGE writes it.
fn genomic_change(s: &str) -> Result<Option<Hit>, String> {
    // Python's `$` also matches before a final newline.
    let s = s.strip_suffix('\n').unwrap_or(s);
    if s.is_empty() || s.chars().any(is_space) {
        return Ok(None);
    }
    // The contig is the longest prefix the rest still matches after.
    for (i, _) in s.rmatch_indices(":g.") {
        if i == 0 {
            continue;
        }
        let rest = &s[i + 3..];
        let digits = rest.find(|c: char| !is_decimal(c)).unwrap_or(rest.len());
        let (num, tail) = rest.split_at(digits);
        let Some((r, a)) = tail.split_once('>') else {
            continue;
        };
        let letters = |x: &str| !x.is_empty() && x.bytes().all(|b| b.is_ascii_alphabetic());
        if !num.is_empty() && letters(r) && letters(a) {
            return Ok(Some(Hit {
                chrom: s[..i].to_owned(),
                pos: int64(num)?,
                reference: r.to_owned(),
                alt: a.to_owned(),
            }));
        }
    }
    Ok(None)
}

/// The variant of a CPSR row, from VAR_ID (`chrom_pos_ref_alt`) or else
/// GENOMIC_CHANGE; `None` when neither reads.
pub fn cpsr_variant(var_id: Option<&str>, change: Option<&str>) -> Result<Option<Hit>, String> {
    let f: Vec<&str> = var_id.unwrap_or_default().split('_').collect();
    let n = f.len();
    if n >= 4 && isdigit(f[n - 3]) {
        return Ok(Some(Hit {
            chrom: f[..n - 3].join("_"),
            pos: int64(f[n - 3])?,
            reference: f[n - 2].to_owned(),
            alt: f[n - 1].to_owned(),
        }));
    }
    genomic_change(change.unwrap_or_default())
}

/// The P/LP rows of a CPSR classification table (final class or CPSR's own) as hits.
pub fn read_cpsr(path: &str) -> Result<Vec<Hit>, String> {
    let t = Table::open(path).map_err(|e| e.to_string())?;
    let mut out = Vec::new();
    for row in t {
        let row = row.map_err(|e| e.to_string())?;
        let plp = [
            "CLASSIFICATION",
            "FINAL_CLASSIFICATION",
            "CPSR_CLASSIFICATION",
        ]
        .iter()
        .any(|c| PLP.contains(&norm_class(row.get(c)).as_str()));
        if let Some(v) = cpsr_variant(row.get("VAR_ID"), row.get("GENOMIC_CHANGE"))?
            && plp
        {
            out.push(v);
        }
    }
    Ok(out)
}

/// A contig name without `chr`, so CPSR's names meet the catalog's.
fn bare(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

type ConstructKey = (String, i64, String, String);

/// The construct catalog: `GENE name` by (contig without chr, pos, REF, ALT).
pub fn read_construct(path: &str) -> Result<HashMap<ConstructKey, String>, String> {
    let mut t: Table = Options::default()
        .skip(|l| l.starts_with('#'))
        .open(path)
        .map_err(|e| e.to_string())?;
    let cols = ["chrom", "pos", "ref", "alt", "gene", "name"];
    let header: Vec<String> = t
        .header()
        .map_err(|e| e.to_string())?
        .map(<[String]>::to_vec)
        .unwrap_or_default();
    let mut out = HashMap::new();
    for row in t {
        let row = row.map_err(|e| e.to_string())?;
        if let Some(miss) = cols.iter().find(|c| !header.iter().any(|x| x == *c)) {
            return Err(format!("{path} has no column {miss}"));
        }
        let pos = row.get("pos").unwrap_or_default();
        if !isdigit(pos) {
            return Err(format!("{path}: position '{pos}' is not a number"));
        }
        out.insert(
            (
                bare(cell(&row, "chrom")?).to_owned(),
                int64(pos)?,
                upper(cell(&row, "ref")?),
                upper(cell(&row, "alt")?),
            ),
            format!("{} {}", cell(&row, "gene")?, cell(&row, "name")?),
        );
    }
    if out.is_empty() {
        return Err(format!("no rows in {path}"));
    }
    Ok(out)
}

/// One output row of the check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckRow {
    pub hit: Hit,
    /// The counts, `None` when the hit was not read.
    pub counts: Option<Counts>,
    pub check: &'static str,
    pub note: String,
    pub construct: String,
}

impl CheckRow {
    fn cells(&self) -> [String; 11] {
        let c = |f: fn(&Counts) -> u64| self.counts.map_or(String::new(), |n| f(&n).to_string());
        [
            self.hit.chrom.clone(),
            self.hit.pos.to_string(),
            self.hit.reference.clone(),
            self.hit.alt.clone(),
            c(|n| n.alt),
            c(|n| n.alt_clipped),
            c(|n| n.reference),
            c(|n| n.ref_clipped),
            self.check.to_owned(),
            self.note.clone(),
            self.construct.clone(),
        ]
    }

    /// The `KNOWN CONSTRUCT SITE` line of the row.
    pub fn construct_line(&self) -> String {
        let r = self.cells();
        let place = format!("{} at {}:{} {}>{}", r[10], r[0], r[1], r[2], r[3]);
        let counts = format!(
            "{} of {} alt reads and {} of {} ref reads clipped at an exon edge",
            r[5], r[4], r[7], r[6]
        );
        if self.check == "not_checked" {
            let tail = if self.counts.is_some() {
                format!(", {counts}")
            } else {
                String::new()
            };
            return format!(
                "  KNOWN CONSTRUCT SITE {place}: not checked ({}){tail}",
                r[9]
            );
        }
        format!("  KNOWN CONSTRUCT SITE {place}: {}, {counts}", r[8])
    }
}

/// The BAM's name for `chrom`: itself, or with `chr` added or removed.
fn bam_contig(chrom: &str, has: impl Fn(&str) -> bool) -> Option<String> {
    [
        Some(chrom.to_owned()),
        Some(format!("chr{chrom}")),
        chrom.strip_prefix("chr").map(str::to_owned),
    ]
    .into_iter()
    .flatten()
    .find(|c| !c.is_empty() && has(c))
}

/// The rows of the check: `skip` writes each hit as not checked with that
/// reason; otherwise each hit's reads are counted from `bam` against `exons`.
pub fn check(
    hits: Vec<Hit>,
    reads: Option<(&str, &str)>,
    skip: Option<&str>,
    construct: Option<&HashMap<ConstructKey, String>>,
) -> Result<Vec<CheckRow>, String> {
    let mut opened = None;
    if skip.is_none()
        && let Some((bam, exons)) = reads
    {
        let ex = Exons::read(exons)?;
        let reader = bam::IndexedReader::from_path(bam).map_err(|e| format!("{bam}: {e}"))?;
        opened = Some((ex, reader));
    }
    let mut out = Vec::with_capacity(hits.len());
    for hit in hits {
        let site_key = (
            bare(&hit.chrom).to_owned(),
            hit.pos,
            upper(&hit.reference),
            upper(&hit.alt),
        );
        let named = construct
            .and_then(|c| c.get(&site_key))
            .cloned()
            .unwrap_or_default();
        let row = |counts, check, note: &str| CheckRow {
            hit: hit.clone(),
            counts,
            check,
            note: note.to_owned(),
            construct: named.clone(),
        };
        let r = match (&mut opened, skip) {
            (_, Some(_)) | (None, None) => row(None, "not_checked", skip.unwrap_or_default()),
            (Some((ex, reader)), None) => {
                if complex_allele(&hit.reference, &hit.alt) {
                    row(None, "not_checked", "not a SNV, MNV or simple indel")
                } else if let Some(contig) =
                    bam_contig(&hit.chrom, |c| reader.header().tid(c.as_bytes()).is_some())
                {
                    let site = Site::new(hit.pos - 1, &hit.reference, &hit.alt);
                    let (lo, hi) = site.window();
                    if lo < 0 {
                        return Err(format!("start out of range ({lo})"));
                    }
                    reader
                        .fetch(FetchDefinition::RegionString(contig.as_bytes(), lo, hi))
                        .map_err(|e| format!("{contig}:{lo}-{hi}: {e}"))?;
                    let reads = reader.records().map(|r| r.map_err(|e| e.to_string()));
                    let n = count_reads(reads, ex, &contig, &site)?;
                    let (check, note) = verdict(n);
                    row(Some(n), check, &note)
                } else {
                    row(
                        None,
                        "not_checked",
                        &format!("the BAM has no contig {}", hit.chrom),
                    )
                }
            }
        };
        out.push(r);
    }
    Ok(out)
}

fn write_table(
    path: &str,
    first: Option<&str>,
    columns: &[&str],
    rows: impl IntoIterator<Item = Vec<String>>,
) -> Result<(), String> {
    let fail = |e: std::io::Error| format!("cannot write {path}: {e}");
    let mut f = AtomicFile::create(path).map_err(fail)?;
    if let Some(line) = first {
        writeln!(f, "{line}").map_err(fail)?;
    }
    writeln!(f, "{}", columns.join("\t")).map_err(fail)?;
    for r in rows {
        writeln!(f, "{}", r.join("\t")).map_err(fail)?;
    }
    f.commit().map_err(fail)
}

fn refuse(err: &mut dyn Write, usage: &str, e: &str) -> u8 {
    let _ = writeln!(err, "{usage}\nog: error: {e}");
    2
}

fn fail(err: &mut dyn Write, e: &str) -> u8 {
    let _ = writeln!(err, "ERROR: {e}");
    2
}

/// `og cdna check ARGS`: the exit code; a refused command line is exit 2.
pub fn check_cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    const USAGE: &str = "usage: og cdna check (--hits HITS | --cpsr CPSR) --out OUT [--bam BAM] \
                         [--exons EXONS] [--skip SKIP] [--construct CONSTRUCT]";
    let a = match Args::parse(
        argv,
        &[
            "--hits",
            "--cpsr",
            "--out",
            "--bam",
            "--exons",
            "--skip",
            "--construct",
        ],
        &[],
    ) {
        Ok(a) => a,
        Err(e) => return refuse(err, USAGE, &e),
    };
    let rows_from = match (a.get("--hits"), a.get("--cpsr")) {
        (Some(_), Some(_)) => {
            return refuse(
                err,
                USAGE,
                "argument --cpsr: not allowed with argument --hits",
            );
        }
        (None, None) => {
            return refuse(err, USAGE, "one of the arguments --hits --cpsr is required");
        }
        (Some(h), None) => Ok(h),
        (None, Some(c)) => Err(c),
    };
    let Some(out_path) = a.get("--out") else {
        return refuse(err, USAGE, "the following arguments are required: --out");
    };
    let skip = a.get("--skip").filter(|s| !s.is_empty());
    let reads = a
        .get("--bam")
        .filter(|s| !s.is_empty())
        .zip(a.get("--exons").filter(|s| !s.is_empty()));
    if skip.is_none() && reads.is_none() {
        return refuse(err, USAGE, "give --bam and --exons, or --skip REASON");
    }
    let run = || -> Result<Vec<CheckRow>, String> {
        let hits = match rows_from {
            Ok(h) => read_hits(h)?,
            Err(c) => read_cpsr(c)?,
        };
        let construct = a
            .get("--construct")
            .filter(|s| !s.is_empty())
            .map(read_construct)
            .transpose()?;
        let rows = check(hits, reads, skip, construct.as_ref())?;
        write_table(
            out_path,
            None,
            &CHECK_COLUMNS,
            rows.iter().map(|r| r.cells().to_vec()),
        )?;
        Ok(rows)
    };
    let rows = match run() {
        Ok(r) => r,
        Err(e) => return fail(err, &e),
    };
    let flagged: Vec<&CheckRow> = rows.iter().filter(|r| r.check == "flagged").collect();
    let clear = rows.iter().filter(|r| r.check == "clear").count();
    let mut said = format!(
        "cDNA read check: {} flagged, {clear} clear, {} not checked\n",
        flagged.len(),
        rows.len() - flagged.len() - clear
    );
    for r in flagged {
        let c = r.cells();
        said += &format!(
            "  FLAGGED {}:{} {}>{}: {} of {} alt reads and {} of {} ref reads clipped at an exon edge ({FLAG_NOTE})\n",
            c[0], c[1], c[2], c[3], c[5], c[4], c[7], c[6]
        );
    }
    for r in rows.iter().filter(|r| !r.construct.is_empty()) {
        said += &r.construct_line();
        said.push('\n');
    }
    match out.write_all(said.as_bytes()) {
        Ok(()) => 0,
        Err(e) => fail(err, &e.to_string()),
    }
}

// --- the scan -----------------------------------------------------------------------

/// One alignment of a read as the scan pairs them: its segment, strand and
/// the read bases it holds, counted in the reference's orientation.
#[derive(Debug, Clone, Copy)]
struct Piece {
    start: i64,
    end: i64,
    left_clip: i64,
    right_clip: i64,
    reverse: bool,
    query_start: i64,
    query_end: i64,
}

fn piece(start: i64, cigar: &[(u8, i64)], reverse: bool) -> Piece {
    let s = segment(start, cigar.iter().copied());
    let held: i64 = cigar
        .iter()
        .filter(|(c, _)| matches!(*c, op::M | op::INS | op::EQ | op::DIFF))
        .map(|(_, n)| n)
        .sum();
    Piece {
        start: s.start,
        end: s.end,
        left_clip: s.left_clip,
        right_clip: s.right_clip,
        reverse,
        query_start: s.left_clip,
        query_end: s.left_clip + held,
    }
}

fn gap_ok(donor: i64, acceptor: i64) -> bool {
    (MIN_INTRON..=MAX_INTRON).contains(&(acceptor - donor))
}

/// The (exon end, exon start) junctions one read shows on `contig`: from
/// its N gaps, and from two of its pieces (the alignment and its SA pieces
/// on this contig) on one strand, the one clipped at an exon end first
/// along the read and the other starting where it stops, give or take
/// [`EDGE_SLACK`] read bases.
fn junctions(
    ex: &Exons,
    contig: &str,
    r: &impl ReadView,
    out: &mut Vec<(i64, i64)>,
) -> Result<(), String> {
    let (ends, starts) = (ex.ends(contig), ex.starts(contig));
    if ends.is_empty() || starts.is_empty() {
        return Ok(());
    }
    let cigar: Vec<(u8, i64)> = ops(r.raw_cigar()).collect();
    let mut rp = r.pos();
    for &(code, n) in &cigar {
        if code == op::SKIP
            && let (Some(d), Some(a)) = (nearest(ends, rp), nearest(starts, rp + n))
            && gap_ok(d, a)
        {
            out.push((d, a));
        }
        if crate::reads::on_reference(code) {
            rp += n;
        }
    }
    let mut pieces = vec![piece(r.pos(), &cigar, r.has(flag::REVERSE))];
    for p in sa_pieces(r.sa().unwrap_or_default()) {
        let p = p?;
        if p.contig == contig {
            pieces.push(piece(p.pos0, &p.cigar, p.reverse));
        }
    }
    for (i, p) in pieces.iter().enumerate() {
        if p.right_clip < MIN_CLIP {
            continue;
        }
        let Some(d) = nearest(ends, p.end) else {
            continue;
        };
        for (j, q) in pieces.iter().enumerate() {
            if i == j
                || q.left_clip < MIN_CLIP
                || q.reverse != p.reverse
                || !(p.query_start < q.query_start
                    && (q.query_start - p.query_end).abs() <= EDGE_SLACK)
            {
                continue;
            }
            if let Some(a) = nearest(starts, q.start)
                && gap_ok(d, a)
            {
                out.push((d, a));
            }
        }
    }
    Ok(())
}

/// A read as a junction counts it: its name and whether it is the first mate.
type ReadKey = (Box<[u8]>, bool);
/// The reads of each (exon end, exon start) junction of one contig.
pub type Junctions = FxHashMap<(i64, i64), FxHashSet<ReadKey>>;

/// Whether a record can show a junction: counted, and with an SA tag or an
/// N gap. Read from the flag, CIGAR and aux bytes, before any decoding.
fn may_show(r: &impl ReadView) -> bool {
    Filter::CDNA.keeps(r)
        && (r.sa().is_some_and(|s| !s.is_empty()) || ops(r.raw_cigar()).any(|(c, _)| c == op::SKIP))
}

/// Adds the junctions of one read to `found`.
pub fn scan_read(
    ex: &Exons,
    contig: &str,
    r: &impl ReadView,
    found: &mut Junctions,
    buf: &mut Vec<(i64, i64)>,
) -> Result<(), String> {
    if !may_show(r) {
        return Ok(());
    }
    buf.clear();
    junctions(ex, contig, r, buf)?;
    for &j in buf.iter() {
        found
            .entry(j)
            .or_default()
            .insert((r.qname().into(), r.has(flag::READ1)));
    }
    Ok(())
}

/// The junctions of one contig, over its merged exon regions. A record that
/// overlaps two regions is read once, in the first: in a later region it
/// starts before the previous region's end.
fn scan_contig(bam: &str, tid: u32, contig: &str, ex: &Exons) -> Result<Junctions, String> {
    let mut reader = bam::IndexedReader::from_path(bam).map_err(|e| format!("{bam}: {e}"))?;
    let mut found = Junctions::default();
    let mut rec = Record::new();
    let mut buf = Vec::new();
    let mut prev_end = i64::MIN;
    for &(a, b) in ex.regions.get(contig).map_or(&[][..], Vec::as_slice) {
        reader
            .fetch(FetchDefinition::Region(
                i32::try_from(tid).map_err(|e| e.to_string())?,
                a,
                b,
            ))
            .map_err(|e| format!("{contig}:{a}-{b}: {e}"))?;
        while let Some(r) = reader.read(&mut rec) {
            r.map_err(|e| format!("{bam}: {e}"))?;
            if rec.pos() < prev_end {
                continue;
            }
            scan_read(ex, contig, &rec, &mut found, &mut buf)?;
        }
        prev_end = b;
    }
    Ok(found)
}

/// One listed locus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub contig: String,
    /// The first junction's exon end, 1-based inclusive (the 0-based end).
    pub start: i64,
    /// The last junction's exon start, 1-based.
    pub end: i64,
    /// Each junction (0-based exon end, 0-based exon start) with its reads.
    pub junctions: Vec<((i64, i64), usize)>,
    /// The locus's distinct reads.
    pub reads: usize,
}

impl Locus {
    fn cells(&self) -> Vec<String> {
        let list: Vec<String> = self
            .junctions
            .iter()
            .map(|&((d, a), n)| format!("{d}-{}:{n}", a + 1))
            .collect();
        vec![
            self.contig.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.junctions.len().to_string(),
            self.reads.to_string(),
            list.join(","),
        ]
    }
}

/// The listed loci of one contig, by position: junctions with
/// [`MIN_READS`] reads, grouped while each donor lies within [`LOCUS_GAP`]
/// of the group's furthest acceptor, kept with [`MIN_JUNCTIONS`] or more.
pub fn loci(contig: &str, found: &Junctions) -> Vec<Locus> {
    let mut kept: Vec<(i64, i64)> = found
        .iter()
        .filter(|(_, k)| k.len() >= MIN_READS)
        .map(|(&j, _)| j)
        .collect();
    kept.sort_unstable();
    let mut groups: Vec<(Vec<(i64, i64)>, i64)> = Vec::new();
    for (d, a) in kept {
        match groups.last_mut() {
            Some((js, far)) if d <= *far + LOCUS_GAP => {
                js.push((d, a));
                *far = (*far).max(a);
            }
            _ => groups.push((vec![(d, a)], a)),
        }
    }
    groups
        .into_iter()
        .filter(|(js, _)| js.len() >= MIN_JUNCTIONS)
        .map(|(js, _)| {
            let reads: FxHashSet<&ReadKey> = js.iter().flat_map(|j| &found[j]).collect();
            Locus {
                contig: contig.to_owned(),
                start: js.iter().map(|j| j.0).min().unwrap_or_default(),
                end: js.iter().map(|j| j.1).max().unwrap_or_default() + 1,
                junctions: js.iter().map(|j| (*j, found[j].len())).collect(),
                reads: reads.len(),
            }
        })
        .collect()
}

/// The scan of a BAM: the loci in the BAM's contig order, the exon regions
/// read and the contigs read. Contigs run on up to `workers` threads.
pub fn scan(bam: &str, ex: &Exons, workers: usize) -> Result<(Vec<Locus>, usize, usize), String> {
    let names: Vec<String> = {
        let r = bam::IndexedReader::from_path(bam).map_err(|e| format!("{bam}: {e}"))?;
        r.header()
            .target_names()
            .iter()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .collect()
    };
    let shared: Vec<(u32, &String)> = (0u32..)
        .zip(&names)
        .filter(|(_, n)| ex.regions.contains_key(n.as_str()))
        .collect();
    if shared.is_empty() {
        return Err(format!(
            "the BAM and the exons have no contig name in common ({bam})"
        ));
    }
    let n_regions = shared
        .iter()
        .map(|(_, n)| ex.regions[n.as_str()].len())
        .sum();
    let next = AtomicUsize::new(0);
    type Slot = Option<Result<Vec<Locus>, String>>;
    let results: Mutex<Vec<Slot>> = Mutex::new(vec![None; shared.len()]);
    thread::scope(|s| {
        for _ in 0..workers.clamp(1, shared.len()) {
            s.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&(tid, name)) = shared.get(i) else {
                        break;
                    };
                    let got = scan_contig(bam, tid, name, ex).map(|f| loci(name, &f));
                    let failed = got.is_err();
                    results.lock().expect("no worker panics")[i] = Some(got);
                    if failed {
                        break;
                    }
                }
            });
        }
    });
    let mut out = Vec::new();
    // A contig left unread lies after a failed one, whose error comes first.
    for r in results
        .into_inner()
        .expect("no worker panics")
        .into_iter()
        .flatten()
    {
        out.extend(r?);
    }
    Ok((out, n_regions, shared.len()))
}

/// `og cdna scan ARGS`: the exit code; a refused command line is exit 2.
pub fn scan_cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    const USAGE: &str = "usage: og cdna scan --out OUT [--bam BAM] [--exons EXONS] [--skip SKIP]";
    let a = match Args::parse(argv, &["--out", "--bam", "--exons", "--skip"], &[]) {
        Ok(a) => a,
        Err(e) => return refuse(err, USAGE, &e),
    };
    let Some(path) = a.get("--out") else {
        return refuse(err, USAGE, "the following arguments are required: --out");
    };
    let write = |comment: &str, rows: &[Locus]| {
        write_table(
            path,
            Some(&format!("# cdna_scan: {comment}")),
            &SCAN_COLUMNS,
            rows.iter().map(Locus::cells),
        )
    };
    if let Some(skip) = a.get("--skip").filter(|s| !s.is_empty()) {
        if let Err(e) = write(&format!("not run: {skip}"), &[]) {
            return fail(err, &e);
        }
        let _ = writeln!(out, "cDNA scan: not run ({skip})");
        return 0;
    }
    let (Some(bam), Some(exons)) = (
        a.get("--bam").filter(|s| !s.is_empty()),
        a.get("--exons").filter(|s| !s.is_empty()),
    ) else {
        return refuse(err, USAGE, "give --bam and --exons, or --skip REASON");
    };
    let run = || -> Result<(Vec<Locus>, usize), String> {
        let ex = Exons::read(exons)?;
        let (rows, n_regions, n_contigs) = scan(bam, &ex, og_io::bgzf::workers().get())?;
        write(
            &format!(
                "read {n_regions} exon regions on {n_contigs} contigs; {} loci listed",
                rows.len()
            ),
            &rows,
        )?;
        Ok((rows, n_regions))
    };
    let (rows, n_regions) = match run() {
        Ok(r) => r,
        Err(e) => return fail(err, &e),
    };
    let mut said = format!(
        "cDNA scan: {} loci with reads across {MIN_JUNCTIONS} or more exon junctions ({n_regions} exon regions read)\n",
        rows.len()
    );
    for r in &rows {
        said += &format!(
            "  SPLICED READS {}:{}-{}: {} exon junctions, {} reads \
             (cDNA or plasmid contamination, or a processed retrocopy: check the reads)\n",
            r.contig,
            r.start,
            r.end,
            r.junctions.len(),
            r.reads
        );
    }
    match out.write_all(said.as_bytes()) {
        Ok(()) => 0,
        Err(e) => fail(err, &e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads::SamRead;
    use crate::reads::tests::{genome, read};

    fn exons_of(spans: &[(i64, i64)], name: &str) -> Exons {
        static N: AtomicUsize = AtomicUsize::new(0);
        let k = N.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("og-cdna-{name}-{}-{k}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("exons.bed");
        let text: String = spans
            .iter()
            .map(|(a, b)| format!("chrT\t{a}\t{b}\n"))
            .collect();
        std::fs::write(&p, text).unwrap();
        let ex = Exons::read(p.to_str().unwrap()).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        ex
    }

    fn sam(name: &str, pos0: i64, cigar: &str, flags: u16, sa: &str, seq: &[u8]) -> SamRead {
        let mut r = read(name, pos0, cigar, seq);
        r.flags = flags;
        r.sa = (!sa.is_empty()).then(|| sa.to_owned());
        r
    }

    // --- the check, on the clip test's fixture -----------------------------------
    const EXON_A: (i64, i64) = (1000, 1200);
    const EXON_B: (i64, i64) = (3000, 3200);
    const POS0: usize = 1150;

    fn alt_base(g: &[u8]) -> u8 {
        b"ACGT"[(b"ACGT".iter().position(|&c| c == g[POS0]).unwrap() + 1) % 4]
    }

    fn with_alt(g: &[u8], seq: &[u8], start: usize) -> Vec<u8> {
        let mut s = seq.to_vec();
        s[POS0 - start] = alt_base(g);
        s
    }

    /// The reads of the three fixtures: spliced, genomic and both.
    fn fixture_reads(kind: &str) -> Vec<SamRead> {
        let g = genome(10_000);
        let mut out = Vec::new();
        let b0 = EXON_B.0 as usize;
        for i in 0..6 {
            let start = 1112 + 2 * i;
            if kind == "both" {
                let aligned = EXON_A.1 as usize - start;
                let seq = [&g[start..EXON_A.1 as usize], &g[b0..b0 + 100 - aligned]].concat();
                out.push(sam(
                    &format!("ref{i}"),
                    start as i64,
                    &format!("{aligned}M{}S", 100 - aligned),
                    0,
                    &format!("chrT,{},+,{aligned}S{}M,60,0;", b0 + 1, 100 - aligned),
                    &seq,
                ));
            } else {
                out.push(sam(
                    &format!("ref{i}"),
                    start as i64,
                    "100M",
                    0,
                    "",
                    &g[start..start + 100],
                ));
            }
        }
        for i in 0..8 {
            let start = 1110 + 2 * i;
            if kind == "genomic" {
                let seq = with_alt(&g, &g[start..start + 100], start);
                out.push(sam(&format!("alt{i}"), start as i64, "100M", 0, "", &seq));
            } else {
                let aligned = EXON_A.1 as usize - start;
                let seq = [
                    with_alt(&g, &g[start..EXON_A.1 as usize], start).as_slice(),
                    &g[b0..b0 + 100 - aligned],
                ]
                .concat();
                out.push(sam(
                    &format!("alt{i}"),
                    start as i64,
                    &format!("{aligned}M{}S", 100 - aligned),
                    0,
                    &format!("chrT,{},+,{aligned}S{}M,60,0;", b0 + 1, 100 - aligned),
                    &seq,
                ));
            }
        }
        if kind != "genomic" {
            out.push(sam(
                "suppl",
                1140,
                "60M40H",
                flag::SUPPLEMENTARY,
                &format!("chrT,{},+,60S40M,60,0;", b0 + 1),
                &with_alt(&g, &g[1140..EXON_A.1 as usize], 1140),
            ));
        }
        out
    }

    fn site() -> Site {
        let g = genome(10_000);
        Site::new(
            POS0 as i64,
            &(g[POS0] as char).to_string(),
            &(alt_base(&g) as char).to_string(),
        )
    }

    fn counts(reads: Vec<SamRead>) -> Counts {
        let ex = exons_of(&[EXON_A, EXON_B], "count");
        count_reads(reads.into_iter().map(Ok), &ex, "chrT", &site()).unwrap()
    }

    fn c(alt: u64, alt_clipped: u64, reference: u64, ref_clipped: u64) -> Counts {
        Counts {
            alt,
            alt_clipped,
            reference,
            ref_clipped,
        }
    }

    #[test]
    fn exon_edges_within_the_slack_on_the_same_side() {
        let ex = exons_of(&[EXON_A, EXON_B], "edges");
        let end = EXON_A.1;
        let right = |e: i64, clip: i64| {
            edge_clipped(&ex, "chrT", e - 100, &[(op::M, 100), (op::SOFT, clip)])
        };
        assert!(right(end, 20));
        assert!(right(end + EDGE_SLACK, 20));
        assert!(!right(end + EDGE_SLACK + 1, 20));
        // The slack is 3 bases, as the Python had it.
        assert!(right(end - 3, 20) && right(end + 3, 20));
        assert!(!right(end - 4, 20) && !right(end + 4, 20));
        assert!(!right(end, MIN_CLIP - 1));
        // A left clip at an exon end (the wrong side) does not count.
        assert!(!edge_clipped(
            &ex,
            "chrT",
            end,
            &[(op::SOFT, 20), (op::M, 100)]
        ));
        assert!(edge_clipped(
            &ex,
            "chrT",
            EXON_B.0,
            &[(op::SOFT, 20), (op::M, 100)]
        ));
        assert!(!edge_clipped(
            &ex,
            "chrQ",
            end - 100,
            &[(op::M, 100), (op::SOFT, 20)]
        ));
        // An SA piece clipped at an exon start counts.
        let r = sam(
            "r",
            1100,
            "90M5S",
            0,
            &format!("chrT,{},+,80S20M,60,0;", EXON_B.0 + 1),
            &[b'A'; 95],
        );
        assert!(read_clipped(&ex, "chrT", &r).unwrap());
        assert_eq!(
            nearest(&[10, 14], 12),
            Some(10),
            "a tie goes to the lower edge"
        );
    }

    #[test]
    fn the_spliced_fixture_counts_and_flags() {
        assert_eq!(counts(fixture_reads("spliced")), c(9, 9, 6, 0));
        assert_eq!(counts(fixture_reads("genomic")), c(8, 0, 6, 0));
        assert_eq!(counts(fixture_reads("both")), c(9, 9, 6, 6));
    }

    #[test]
    fn a_repeated_name_a_duplicate_and_a_secondary_alignment_do_not_count() {
        let g = genome(10_000);
        let seq = with_alt(&g, &g[1110..1210], 1110);
        let mut reads = fixture_reads("spliced");
        reads.push(sam("alt0", 1110, "100M", 0, "", &seq));
        reads.push(sam("dup", 1110, "100M", flag::DUPLICATE, "", &seq));
        reads.push(sam("sec", 1110, "100M", flag::SECONDARY, "", &seq));
        assert_eq!(counts(reads), c(9, 9, 6, 0));
    }

    #[test]
    fn the_verdicts() {
        assert_eq!(verdict(c(9, 9, 6, 0)), ("flagged", FLAG_NOTE.to_owned()));
        assert_eq!(verdict(c(9, 9, 6, 6)).0, "clear");
        assert_eq!(verdict(c(9, 2, 6, 0)).0, "clear");
        assert_eq!(verdict(c(9, 9, 0, 0)).0, "flagged");
        assert_eq!(
            verdict(c(1, 1, 6, 0)),
            (
                "not_checked",
                "fewer than 2 reads carry the alt allele".to_owned()
            )
        );
        // 0.2 * 6 is 1.2000000000000002 in floating point: 1 of 6 still flags.
        assert_eq!(verdict(c(4, 2, 6, 1)).0, "flagged");
        assert_eq!(verdict(c(4, 2, 10, 2)).0, "flagged");
        assert_eq!(verdict(c(4, 2, 10, 3)).0, "clear");
    }

    #[test]
    fn the_construct_lines() {
        let hit = Hit {
            chrom: "chrT".into(),
            pos: 1151,
            reference: "A".into(),
            alt: "C".into(),
        };
        let row = |counts, check, note: &str| CheckRow {
            hit: hit.clone(),
            counts,
            check,
            note: note.into(),
            construct: "GENET X1Y".into(),
        };
        assert_eq!(
            row(Some(c(8, 0, 6, 0)), "clear", "").construct_line(),
            "  KNOWN CONSTRUCT SITE GENET X1Y at chrT:1151 A>C: clear, 0 of 8 alt reads and 0 of 6 ref reads \
             clipped at an exon edge"
        );
        assert_eq!(
            row(
                Some(c(1, 0, 6, 0)),
                "not_checked",
                "fewer than 2 reads carry the alt allele"
            )
            .construct_line(),
            "  KNOWN CONSTRUCT SITE GENET X1Y at chrT:1151 A>C: not checked (fewer than 2 reads carry the alt \
             allele), 0 of 1 alt reads and 0 of 6 ref reads clipped at an exon edge"
        );
        assert_eq!(
            row(None, "not_checked", "no BAM").construct_line(),
            "  KNOWN CONSTRUCT SITE GENET X1Y at chrT:1151 A>C: not checked (no BAM)"
        );
    }

    #[test]
    fn the_shipped_construct_catalog() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../../assets/construct_variants.tsv"
        );
        let cat = read_construct(path).unwrap();
        assert_eq!(cat.len(), 10);
        assert!(cat.keys().all(|k| k.1 > 0
            && k.2.chars().all(|c| c.is_ascii_alphabetic())
            && k.3.chars().all(|c| c.is_ascii_alphabetic())
            && k.2 != k.3));
        let key = ("12".to_owned(), 40_340_400, "G".to_owned(), "A".to_owned());
        assert_eq!(cat.get(&key).map(String::as_str), Some("LRRK2 G2019S"));
    }

    #[test]
    fn cpsr_variants_and_classes() {
        let v = |id: Option<&str>, ch: Option<&str>| {
            cpsr_variant(id, ch)
                .unwrap()
                .map(|h| (h.chrom, h.pos, h.reference, h.alt))
        };
        assert_eq!(
            v(Some("13_32340300_G_A"), None),
            Some(("13".into(), 32_340_300, "G".into(), "A".into()))
        );
        assert_eq!(
            v(Some("HLA_A_5_C_T"), None),
            Some(("HLA_A".into(), 5, "C".into(), "T".into()))
        );
        assert_eq!(
            v(Some(""), Some("2:g.47800000A>G")),
            Some(("2".into(), 47_800_000, "A".into(), "G".into()))
        );
        assert_eq!(
            v(None, Some("chrUn:g.x:g.12AC>T\n")),
            Some(("chrUn:g.x".into(), 12, "AC".into(), "T".into()))
        );
        assert_eq!(v(None, Some("2:g.47800000A>G T")), None);
        assert_eq!(v(None, Some("2:g.A>G")), None);
        assert_eq!(norm_class(Some(" Likely_Pathogenic ")), "likely pathogenic");
        assert_eq!(norm_class(Some("LP")), "likely pathogenic");
        assert_eq!(norm_class(Some("NA")), "");
        assert_eq!(norm_class(None), "");
    }

    // --- the scan, on the scan test's fixture ------------------------------------
    const E1: (i64, i64) = (1000, 1200);
    const E2: (i64, i64) = (3000, 3150);
    const E3: (i64, i64) = (5000, 5200);

    fn across(name: &str, left: (i64, i64), right: (i64, i64), aligned: i64) -> SamRead {
        let clip = 100 - aligned;
        sam(
            name,
            left.1 - aligned,
            &format!("{aligned}M{clip}S"),
            0,
            &format!("chrT,{},+,{aligned}S{clip}M,60,0;", right.0 + 1),
            &[b'A'; 100],
        )
    }

    fn scan_fixture(kind: &str) -> Vec<SamRead> {
        let a = |i: i64| 60 + 5 * i;
        match kind {
            "genomic" => {
                let mut v = Vec::new();
                for e in [E1, E2, E3] {
                    for i in 0..4 {
                        v.push(sam(
                            &format!("g{i}"),
                            e.1 - 50 + 7 * i,
                            "100M",
                            0,
                            "",
                            &[b'A'; 100],
                        ));
                    }
                }
                for e in [E2, E3] {
                    for i in 0..4 {
                        v.push(sam(
                            &format!("s{i}"),
                            e.0 - 50 + 7 * i,
                            "100M",
                            0,
                            "",
                            &[b'A'; 100],
                        ));
                    }
                }
                v
            }
            "one" => (0..4)
                .map(|i| across(&format!("a{i}"), E1, E2, a(i)))
                .collect(),
            "thin" => vec![
                across("a0", E1, E2, 60),
                across("b0", E2, E3, 60),
                across("c0", E1, E3, 60),
            ],
            _ => (0..3)
                .map(|i| across(&format!("a{i}"), E1, E2, a(i)))
                .chain([sam(
                    "n0",
                    E1.1 - 50,
                    &format!("50M{}N50M", E2.0 - E1.1),
                    0,
                    "",
                    &[b'A'; 100],
                )])
                .chain((0..3).map(|i| across(&format!("b{i}"), E2, E3, a(i))))
                .chain((0..2).map(|i| across(&format!("c{i}"), E1, E3, a(i))))
                .collect(),
        }
    }

    fn found(ex: &Exons, reads: &[SamRead]) -> Junctions {
        let mut f = Junctions::default();
        let mut buf = Vec::new();
        for r in reads {
            scan_read(ex, "chrT", r, &mut f, &mut buf).unwrap();
        }
        f
    }

    fn listed(kind: &str, extra: &[SamRead]) -> Vec<Vec<String>> {
        let ex = exons_of(&[E1, E2, E3], &format!("scan-{kind}"));
        let mut reads = scan_fixture(kind);
        reads.extend_from_slice(extra);
        loci("chrT", &found(&ex, &reads))
            .iter()
            .map(Locus::cells)
            .collect()
    }

    fn spliced_row() -> Vec<String> {
        vec![
            "chrT".into(),
            E1.1.to_string(),
            (E3.0 + 1).to_string(),
            "3".into(),
            "9".into(),
            format!(
                "{}-{}:4,{}-{}:2,{}-{}:3",
                E1.1,
                E2.0 + 1,
                E1.1,
                E3.0 + 1,
                E2.1,
                E3.0 + 1
            ),
        ]
    }

    fn junctions_of(sa: &str, flags: u16) -> Vec<(i64, i64)> {
        let ex = exons_of(&[E1, E2, E3], "junctions");
        let r = sam("p", E1.1 - 60, "60M40S", flags, sa, &[b'A'; 100]);
        let mut out = Vec::new();
        junctions(&ex, "chrT", &r, &mut out).unwrap();
        out.sort_unstable();
        out.dedup();
        out
    }

    #[test]
    fn merged_regions_and_the_nearest_edge() {
        let ex = exons_of(&[E1, E2, E3, (1100, 1300), (1300, 1310)], "merge");
        assert_eq!(ex.regions["chrT"], vec![(1000, 1310), E2, E3]);
        assert_eq!(nearest(ex.ends("chrT"), E1.1 + EDGE_SLACK), Some(E1.1));
        assert_eq!(nearest(&[E1.1], E1.1 + EDGE_SLACK + 1), None);
    }

    #[test]
    fn the_scan_lists_spliced_genes_only() {
        assert_eq!(listed("spliced", &[]), vec![spliced_row()]);
        assert!(listed("genomic", &[]).is_empty(), "genomic reads only");
        assert!(listed("one", &[]).is_empty(), "one junction alone");
        assert!(listed("thin", &[]).is_empty(), "one read per junction");
    }

    #[test]
    fn reads_that_do_not_count_twice_or_at_all() {
        let a = across("a0", E1, E2, 60);
        let mut d = a.clone();
        d.qname = b"d".to_vec();
        d.flags = flag::DUPLICATE;
        let mut s = a.clone();
        s.qname = b"s".to_vec();
        s.flags = flag::SECONDARY;
        let sup = sam(
            "a0",
            E2.0,
            "60S40M",
            flag::SUPPLEMENTARY,
            &format!("chrT,{},+,60M40S,60,0;", a.pos + 1),
            &[b'A'; 100],
        );
        assert_eq!(listed("spliced", &[a, d, s, sup]), vec![spliced_row()]);
    }

    #[test]
    fn a_read_with_its_primary_elsewhere_counts_from_its_supplementary_records() {
        let ex = exons_of(&[E1, E2, E3], "vec");
        let v1 = sam(
            "v",
            E1.1 - 50,
            "50S50M50S",
            flag::SUPPLEMENTARY,
            &format!("chrV,1,+,50M100S,60,0;chrT,{},+,100S50M,60,0;", E2.0 + 1),
            &[b'A'; 150],
        );
        let v2 = sam(
            "v",
            E2.0,
            "100S50M",
            flag::SUPPLEMENTARY,
            &format!("chrV,1,+,50M100S,60,0;chrT,{},+,50S50M50S,60,0;", E1.1 - 49),
            &[b'A'; 150],
        );
        let f = found(&ex, &[v1, v2]);
        assert_eq!(f.len(), 1);
        assert_eq!(f[&(E1.1, E2.0)].len(), 1);
    }

    #[test]
    fn pieces_pair_on_one_strand_in_read_order() {
        let j = vec![(E1.1, E2.0)];
        let b = E2.0 + 1;
        assert_eq!(junctions_of(&format!("chrT,{b},+,60S40M,60,0;"), 0), j);
        let s = EDGE_SLACK;
        assert_eq!(
            junctions_of(&format!("chrT,{b},+,{}S{}M,60,0;", 60 - s, 40 + s), 0),
            j
        );
        assert_eq!(junctions_of(&format!("chrT,{b},-,60S40M,60,0;"), 0), vec![]);
        assert_eq!(
            junctions_of(&format!("chrT,{b},+,10S40M50S,60,0;"), 0),
            vec![]
        );
        assert_eq!(
            junctions_of(&format!("chrT,{b},-,60S40M,60,0;"), flag::REVERSE),
            j
        );
        assert_eq!(
            junctions_of(&format!("chrT,{b},+,{}S{}M,60,0;", 60 + s, 40 - s), 0),
            j
        );
        assert_eq!(
            junctions_of(&format!("chrT,{b},+,{}S{}M,60,0;", 61 + s, 39 - s), 0),
            vec![]
        );
    }

    #[test]
    fn a_read_across_two_junctions_counts_on_each_and_once_in_the_locus() {
        let two = sam(
            "n1",
            E1.1 - 50,
            &format!("50M{}N150M{}N50M", E2.0 - E1.1, E3.0 - E2.1),
            0,
            "",
            &[b'A'; 250],
        );
        let rows = listed("spliced", &[two]);
        assert_eq!(rows.len(), 1);
        assert_eq!(
            rows[0][4..],
            [
                "10".to_owned(),
                format!(
                    "{}-{}:5,{}-{}:2,{}-{}:4",
                    E1.1,
                    E2.0 + 1,
                    E1.1,
                    E3.0 + 1,
                    E2.1,
                    E3.0 + 1
                )
            ]
        );
        assert!(gap_ok(1200, 1200 + MIN_INTRON) && !gap_ok(1200, 1199 + MIN_INTRON));
    }

    #[test]
    fn a_locus_needs_two_junctions_close_together() {
        let ex = exons_of(&[E1, E2, E3], "far");
        let key = |n: &str| (n.as_bytes().into(), true);
        let mut f = Junctions::default();
        for (j, names) in [
            ((100, 300), ["a", "b"]),
            ((400, 600), ["b", "c"]),
            ((600 + LOCUS_GAP + 1, 700 + LOCUS_GAP), ["d", "e"]),
        ] {
            f.insert(j, names.iter().map(|n| key(n)).collect());
        }
        let got = loci("chrT", &f);
        assert_eq!(got.len(), 1, "the far junction stands alone");
        assert_eq!((got[0].start, got[0].end, got[0].reads), (100, 601, 3));
        drop(ex);
    }
}
