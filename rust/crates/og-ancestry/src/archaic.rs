//! Step 43, archaic ancestry (`docs/43-archaic.md`): the commands around
//! IBDmix (Chen et al. 2020), as `og archaic <verb>`.
//!
//! IBDmix compares a modern genome with an archaic one site by site and
//! adds up a log-odds score (LOD) that a stretch of the modern genome shares
//! a recent ancestor with the archaic genome. It takes each site's allele
//! frequency from the modern people in its input, and one genome alone
//! gives no frequency, so `genotypes` adds stand-in columns whose genotypes
//! carry the 1000 Genomes allele count of the chosen population.
//!
//! - `table` (setup): one GRCh37 chromosome of the archaic table, a
//!   sites-only VCF that Picard lifts to GRCh38. The 1000 Genomes SNVs and
//!   the three archaic genomes' calls are merged by position in one pass;
//!   each genome is read inside its own mask, and a call outside it (or at
//!   an excluded site) is kept as masked, code 10 plus its ALT count.
//! - `strict` (setup): the lifted rows of one chromosome, once per
//!   position, with every call outside the 1000 Genomes strict mask marked
//!   masked.
//! - `genotypes` (step 43): IBDmix's genotype file of each chromosome, each
//!   archaic genome's mask of masked positions, and the spans table.
//! - `report` (step 43): the segments at LOD 4 and 50 kb or longer, merged
//!   into Neanderthal and Denisovan ones, the shares, the ten longest and
//!   their protein-coding genes.
//!
//! Each command writes what `bin/archaic_report.py` wrote (checked by
//! `tests/parity`). Deliberate changes: every output is written under a
//! temporary name and renamed into place; a gzip input is told by its first
//! bytes, not its name; `strict` reads its rows in position order and stops
//! (exit 2) on a row that goes back, instead of keeping every position of a
//! chromosome in memory; and malformed input is an `ERROR:` line, never a
//! traceback.
//!
//! The table's bytes are part of the archaic data version: a change to what
//! `table` or `strict` writes must change `ARCHAIC_TABLE_FORMAT` in
//! `scripts/setup.sh`, so existing installs rebuild their table.

use super::{ArgError, Args, Fail, Result, Spec, cannot_read, cannot_write, int, repr, split_ws};
use flate2::Compression;
use flate2::write::GzEncoder;
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::fmt::{fixed, general, zero_pad};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

pub const ARCHAICS: [&str; 3] = ["Altai", "Vindija", "Denisova"];
const NEANDERTHALS: [&str; 2] = ["Altai", "Vindija"];
const POPS: [&str; 6] = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS"];
const TABLE_HEADER: [&str; 13] = [
    "chrom", "pos", "ref", "alt", "Altai", "Vindija", "Denisova", "AF_ALL", "AF_AFR", "AF_AMR",
    "AF_EAS", "AF_EUR", "AF_SAS",
];
/// Table code: MASKED plus the ALT count, a call outside its genome's mask.
const MASKED: u8 = 10;
/// Stand-in modern people: 200 alleles, about one 1000 Genomes population.
const PSEUDO: i64 = 100;
/// Chen et al. 2020: LOD 4 and 50 kb.
const MIN_LOD: f64 = 4.0;
const MIN_LENGTH: i64 = 50_000;
const TOP: usize = 10;
const SEGMENT_HEADER: &str = "chrom\tstart\tend\tlength_bp\tancestry\tmatched\tlod\tgenes";
const CAVEAT: &str = "IBDmix finds stretches where your genome matches a Neanderthal or the Denisovan genome more \
closely than the allele frequencies of a modern population explain. A stretch counts when either \
of your two copies carries it, and the share is of the analysed span of the autosomes. Some of what \
it finds in anyone is old ancestry shared with the archaic genomes rather than gene flow from them: \
Chen et al. 2020 find Neanderthal sequence in every population they studied, Africans included. \
The numbers move with the cutoffs (LOD 4, 50 kb) and with the frequencies used. An estimate, not \
a measurement.";

/// Why a command stopped: its exit code (2 for a wrong option value the
/// command checks itself, 1 otherwise) and the message after `ERROR: `.
#[derive(Debug)]
struct Stop(u8, String);

impl From<Fail> for Stop {
    fn from(f: Fail) -> Self {
        Stop(1, f.0)
    }
}

fn die<T>(msg: impl Into<String>) -> std::result::Result<T, Stop> {
    Err(Stop(1, msg.into()))
}

fn usage<T>(msg: impl Into<String>) -> std::result::Result<T, Stop> {
    Err(Stop(2, msg.into()))
}

type Out<T> = std::result::Result<T, Stop>;

/// A text file read line by line as Python's text mode reads it (universal
/// newlines, strict UTF-8; gzip read through), keeping the raw line for
/// messages that quote it.
struct Text {
    path: PathBuf,
    lines: Lines<Box<dyn BufRead>>,
    line: String,
}

impl Text {
    fn open(path: &str) -> Result<Text> {
        let path = PathBuf::from(path);
        let input = og_io::open(&path).map_err(|e| cannot_read(&path, &e))?;
        Ok(Text {
            path,
            lines: Lines::new(input, Utf8::Strict),
            line: String::new(),
        })
    }

    /// The next line with its `\n`, or `None` at the end.
    fn raw(&mut self) -> Result<Option<&str>> {
        let more = self
            .lines
            .read_line(&mut self.line)
            .map_err(|e| cannot_read(&self.path, &e))?;
        Ok(more.then_some(self.line.as_str()))
    }

    /// The next line without its `\n` (`line.rstrip("\n")`).
    fn next(&mut self) -> Result<Option<&str>> {
        Ok(self.raw()?.map(|l| l.strip_suffix('\n').unwrap_or(l)))
    }
}

/// Python's `int(s)`, or the message of the ValueError it raises.
fn need_int(s: &str) -> Out<i64> {
    int(s).ok_or_else(|| {
        Stop(
            1,
            format!("invalid literal for int() with base 10: {}", repr(s)),
        )
    })
}

/// Python's `float(s)`, or the message of the ValueError it raises.
fn need_float(s: &str) -> Out<f64> {
    pycompat::num::parse_float(s)
        .map_err(|_| Stop(1, format!("could not convert string to float: {}", repr(s))))
}

/// `repr(line[:80])`.
fn head80(line: &str) -> String {
    let end = line.char_indices().nth(80).map_or(line.len(), |(i, _)| i);
    repr(&line[..end])
}

/// `os.path.join(a, b)` for a relative `b`.
fn join(a: &str, b: &str) -> String {
    if a.is_empty() || a.ends_with('/') {
        format!("{a}{b}")
    } else {
        format!("{a}/{b}")
    }
}

/// `'0/1'`, `'1|1'`, `'0/0'` to 1, 2, 0; a missing, partial or
/// multi-allelic call to `None`.
fn gt_code(gt: &str) -> Option<u8> {
    let mut parts = gt.split(['/', '|']);
    let (Some(a), Some(b), None) = (parts.next(), parts.next(), parts.next()) else {
        return None;
    };
    let allele = |x: &str| match x {
        "0" => Some(0),
        "1" => Some(1),
        _ => None,
    };
    Some(allele(a)? + allele(b)?)
}

/// Python's `max(a, b)` and `min(a, b)` of floats: the first unless the
/// second compares greater (smaller), so a NaN first stays.
fn py_max(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

fn py_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

// --- table (setup) ----------------------------------------------------------------

/// One mask's intervals on one chromosome, sorted and merged: 0-based,
/// half-open, so `[s, e)` holds the 1-based positions `s + 1` to `e`. The
/// chromosome is matched with and without the `chr` prefix.
fn read_mask(path: &str, chrom: &str) -> Out<Vec<(i64, i64)>> {
    let num = chrom.strip_prefix("chr").unwrap_or(chrom);
    let with_chr = format!("chr{num}");
    let mut iv = Vec::new();
    let mut t = Text::open(path)?;
    while let Some(line) = t.raw()? {
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let c: Vec<&str> = split_ws(line).collect();
        if c.len() < 3 || (c[0] != num && c[0] != with_chr) {
            continue;
        }
        let (s, e) = (need_int(c[1])?, need_int(c[2])?);
        if e <= s {
            continue; // an empty or reversed interval holds no position
        }
        iv.push((s, e));
    }
    iv.sort_unstable();
    let mut out: Vec<(i64, i64)> = Vec::with_capacity(iv.len());
    for (s, e) in iv {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    Ok(out)
}

/// Whether a 1-based position is inside a mask, for positions asked in
/// ascending order: the walk over the intervals never goes back.
struct Mask {
    iv: Vec<(i64, i64)>,
    i: usize,
}

impl Mask {
    fn holds(&mut self, pos: i64) -> bool {
        let p = pos - 1;
        while self.i < self.iv.len() && self.iv[self.i].1 <= p {
            self.i += 1;
        }
        self.i < self.iv.len() && self.iv[self.i].0 <= p
    }
}

/// One record of a table input: a 1000 Genomes SNV (REF, ALT and the six
/// frequencies as written) or an archaic call (REF, ALT, GT).
#[derive(Debug, Clone)]
enum Rec {
    Modern {
        r#ref: String,
        alt: String,
        afs: Vec<String>,
    },
    Archaic {
        r#ref: String,
        alt: String,
        gt: String,
    },
}

/// One input of the merge: its lines, read on demand.
struct Stream {
    text: Text,
    /// 0 for the 1000 Genomes lines, 1 + i for archaic genome i.
    src: usize,
}

impl Stream {
    /// The next record, `(pos, rec)`; 1000 Genomes lines that are not
    /// SNVs are skipped.
    fn next(&mut self) -> Out<Option<(i64, Rec)>> {
        loop {
            let Some(raw) = self.text.raw()? else {
                return Ok(None);
            };
            let line = raw.strip_suffix('\n').unwrap_or(raw);
            let c: Vec<&str> = line.split('\t').collect();
            if self.src == 0 {
                if c.len() < 9 {
                    let (n, head) = (c.len(), head80(raw));
                    return die(format!(
                        "{}: a line has {n} columns, want 9: {head}",
                        self.text.path.display()
                    ));
                }
                if c[1].chars().count() != 1 || c[2].chars().count() != 1 {
                    continue;
                }
                let pos = need_int(c[0])?;
                return Ok(Some((
                    pos,
                    Rec::Modern {
                        r#ref: c[1].to_owned(),
                        alt: c[2].to_owned(),
                        afs: c[3..9].iter().map(|s| (*s).to_owned()).collect(),
                    },
                )));
            }
            if c.len() < 4 {
                let (n, head) = (c.len(), head80(raw));
                return die(format!(
                    "{}: a line has {n} columns, want 4: {head}",
                    self.text.path.display()
                ));
            }
            let pos = need_int(c[0])?;
            return Ok(Some((
                pos,
                Rec::Archaic {
                    r#ref: c[1].to_owned(),
                    alt: c[2].to_owned(),
                    gt: c[3].to_owned(),
                },
            )));
        }
    }
}

/// An archaic call at one position: (REF, ALT, GT).
type Call<'a> = Option<(&'a str, &'a str, &'a str)>;

/// The ALT count of one archaic record at a row with this REF and ALT, or 9.
fn archaic_code(rec: Call<'_>, r: &str, alt: &str) -> u8 {
    let Some((aref, aalt, gt)) = rec else {
        return 9;
    };
    if aref != r {
        return 9;
    }
    let Some(code) = gt_code(gt) else {
        return 9;
    };
    if aalt == "." || aalt == "<*>" || code == 0 {
        return if code == 0 { 0 } else { 9 };
    }
    if aalt == alt { code } else { 9 }
}

/// The row of one position: REF, ALT, the three codes and the six
/// frequencies as written, or `None` when nothing is known there.
fn table_row<'a>(
    modern: Option<(&'a str, &'a str, &'a [String])>,
    arch: [Call<'a>; 3],
) -> Option<(&'a str, &'a str, [u8; 3], Vec<String>)> {
    const ZEROS: [&str; 6] = ["0"; 6];
    let (r, alt, afs): (&str, &str, Vec<&str>) = match modern {
        Some((r, alt, afs)) => (r, alt, afs.iter().map(String::as_str).collect()),
        None => {
            // A site no modern person varies at: the ALT of the first
            // archaic genome with a non-reference call, every frequency 0.
            let (r, alt, _) = arch.iter().flatten().find(|(r, a, gt)| {
                !matches!(gt_code(gt), None | Some(0))
                    && r.chars().count() == 1
                    && a.chars().count() == 1
            })?;
            (r, alt, ZEROS.to_vec())
        }
    };
    let codes = arch.map(|a| archaic_code(a, r, alt));
    if codes.iter().all(|&c| c == 9) {
        return None;
    }
    if modern.is_none() && !codes.iter().any(|&c| c == 1 || c == 2) {
        return None;
    }
    let clean = afs
        .iter()
        .map(|v| {
            let x = pycompat::num::parse_float(v).unwrap_or(0.0);
            general(py_min(py_max(x, 0.0), 1.0), 6)
        })
        .collect();
    Some((r, alt, codes, clean))
}

/// `NAME=FILE` with NAME one of the archaic genomes.
fn named<'a>(spec: &'a str, option: &str) -> Out<(&'a str, &'a str)> {
    let (name, path) = spec.split_once('=').unwrap_or((spec, ""));
    if !ARCHAICS.contains(&name) || path.is_empty() {
        return usage(format!(
            "{option} takes NAME=FILE with NAME one of {}, got {}",
            ARCHAICS.join(", "),
            repr(spec)
        ));
    }
    Ok((name, path))
}

fn cmd_table(a: &Args, err: &mut dyn Write) -> Out<()> {
    let chrom_arg = a.req("--chrom");
    let mut names = Vec::new();
    let mut paths = Vec::new();
    for spec in a.all("--archaic") {
        let (name, path) = named(spec, "--archaic")?;
        names.push(name);
        paths.push(path);
    }
    if names != ARCHAICS {
        return usage(format!(
            "--archaic must name {} in that order, got {}",
            ARCHAICS.join(", "),
            names.join(", ")
        ));
    }
    let mut masks: [Option<Mask>; 3] = [None, None, None];
    for spec in a.all("--mask") {
        let (name, path) = named(spec, "--mask")?;
        let k = ARCHAICS.iter().position(|n| *n == name).expect("checked");
        if masks[k].is_some() {
            return usage(format!("--mask names {name} twice"));
        }
        let iv = read_mask(path, chrom_arg)?;
        if iv.is_empty() {
            return die(format!(
                "{path} holds no interval on chromosome {chrom_arg}: a mask that matches nothing \
                 would make every site unknown"
            ));
        }
        masks[k] = Some(Mask { iv, i: 0 });
    }
    let mut excluded = Vec::new();
    for path in a.all("--exclude") {
        let iv = read_mask(path, chrom_arg)?;
        if iv.is_empty() {
            return die(format!(
                "{path} holds no interval on chromosome {chrom_arg}: an empty exclusion list is a \
                 download that went wrong, not a chromosome without segmental duplications"
            ));
        }
        excluded.push(Mask { iv, i: 0 });
    }
    let chrom = if chrom_arg.starts_with("chr") {
        chrom_arg.to_owned()
    } else {
        format!("chr{chrom_arg}")
    };
    let out_path = PathBuf::from(a.req("--out"));
    let mut out = AtomicFile::create(&out_path).map_err(|e| cannot_write(&out_path, &e))?;
    let mut head = format!("##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n");
    for n in ARCHAICS {
        head += &format!(
            "##INFO=<ID={n},Number=1,Type=Integer,Description=\"ALT alleles in the {n} genome (9: unknown; 10 to 12: 0 to 2 outside its mask)\">\n"
        );
    }
    for p in POPS {
        head += &format!(
            "##INFO=<ID=AF_{p},Number=1,Type=Float,Description=\"1000 Genomes phase 3 ALT frequency ({p})\">\n"
        );
    }
    head += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let w = |out: &mut AtomicFile, s: &str| {
        out.write_all(s.as_bytes())
            .map_err(|e| Stop::from(cannot_write(&out_path, &e)))
    };
    w(&mut out, &head)?;

    // The k-way merge of heapq.merge: the head of each input on a heap
    // keyed (position, input), each input advanced only when its head is
    // taken, so the order is the same as the Python's on any input.
    let mut streams = vec![Stream {
        text: Text::open(a.req("--modern"))?,
        src: 0,
    }];
    for (i, p) in paths.iter().enumerate() {
        streams.push(Stream {
            text: Text::open(p)?,
            src: 1 + i,
        });
    }
    let mut heads: Vec<Option<Rec>> = vec![None; streams.len()];
    let mut heap = BinaryHeap::new();
    for (i, s) in streams.iter_mut().enumerate() {
        if let Some((pos, rec)) = s.next()? {
            heads[i] = Some(rec);
            heap.push(Reverse((pos, i)));
        }
    }
    let mut rows: u64 = 0;
    let mut cur: Option<i64> = None;
    let mut modern: Option<Rec> = None;
    let mut arch: [Option<Rec>; 3] = [None, None, None];
    let mut last = -1i64;
    let mut flush = |cur: i64,
                     modern: &Option<Rec>,
                     arch: &[Option<Rec>; 3],
                     out: &mut AtomicFile|
     -> Out<()> {
        // Every mask is walked for every position, before any early
        // return, because the walk only goes forward.
        let mut excl = false;
        for m in excluded.iter_mut() {
            excl |= m.holds(cur);
        }
        let mut outside = [false; 3];
        for (k, m) in masks.iter_mut().enumerate() {
            outside[k] = excl || m.as_mut().is_some_and(|m| !m.holds(cur));
        }
        let modern = modern.as_ref().map(|m| match m {
            Rec::Modern { r#ref, alt, afs } => (r#ref.as_str(), alt.as_str(), afs.as_slice()),
            Rec::Archaic { .. } => unreachable!("only the modern input gives Modern"),
        });
        let calls = [0, 1, 2].map(|k| match &arch[k] {
            Some(Rec::Archaic { r#ref, alt, gt }) => {
                Some((r#ref.as_str(), alt.as_str(), gt.as_str()))
            }
            _ => None,
        });
        let Some((r, alt, codes, afs)) = table_row(modern, calls) else {
            return Ok(());
        };
        // A call outside its genome's mask, or at an excluded site, is kept
        // but marked masked: step 43 hands IBDmix the call and a mask of
        // those positions, as Chen et al. do.
        let mut info = String::new();
        for k in 0..3 {
            let c = codes[k];
            let c = if outside[k] && c != 9 { c + MASKED } else { c };
            info += &format!("{}={c};", ARCHAICS[k]);
        }
        for (k, p) in POPS.iter().enumerate() {
            info += &format!("AF_{p}={}", afs[k]);
            if k + 1 < POPS.len() {
                info.push(';');
            }
        }
        out.write_all(format!("{chrom}\t{cur}\t.\t{r}\t{alt}\t.\t.\t{info}\n").as_bytes())
            .map_err(|e| Stop::from(cannot_write(&out_path, &e)))?;
        rows += 1;
        Ok(())
    };
    while let Some(Reverse((pos, i))) = heap.pop() {
        let rec = heads[i].take().expect("a head for each heap entry");
        if pos < last {
            return die(format!("input not sorted by position at {pos}"));
        }
        if cur != Some(pos) {
            if let Some(c) = cur {
                flush(c, &modern, &arch, &mut out)?;
            }
            cur = Some(pos);
            modern = None;
            arch = [None, None, None];
        }
        last = pos;
        match rec {
            Rec::Modern { .. } => {
                // A second modern record at one position (a split
                // multi-allelic site) makes the row ambiguous: keep neither.
                modern = Some(match modern {
                    None => rec,
                    Some(_) => Rec::Modern {
                        r#ref: "N".into(),
                        alt: "N".into(),
                        afs: vec!["0".into(); POPS.len()],
                    },
                });
            }
            Rec::Archaic { .. } => arch[streams[i].src - 1] = Some(rec),
        }
        // As heapq.merge does, the input whose record was taken is read on
        // only now.
        if let Some((p, r)) = streams[i].next()? {
            heads[i] = Some(r);
            heap.push(Reverse((p, i)));
        }
    }
    if let Some(c) = cur {
        flush(c, &modern, &arch, &mut out)?;
    }
    out.commit().map_err(|e| cannot_write(&out_path, &e))?;
    let _ = writeln!(err, "{chrom}: {rows} rows");
    Ok(())
}

// --- strict (setup) ---------------------------------------------------------------

/// The lifted GRCh38 rows of one chromosome (no header, the 13 table
/// columns), once per position, with the calls of every row outside the
/// 1000 Genomes strict mask marked masked (MASKED plus the ALT count; 9
/// stays 9). The strict mask is the same for the three genomes, so all
/// three are marked.
fn cmd_strict(a: &Args, err: &mut dyn Write) -> Out<()> {
    let chrom = a.req("--chrom");
    let strict = a.req("--strict");
    let iv = read_mask(strict, chrom)?;
    if iv.is_empty() {
        return die(format!(
            "{strict} holds no interval on chromosome {chrom}: a mask that matches nothing \
             would mark every site masked"
        ));
    }
    let rows_path = a.req("--rows");
    let mut t = Text::open(rows_path)?;
    let out_path = PathBuf::from(a.req("--out"));
    let mut out = AtomicFile::create(&out_path).map_err(|e| cannot_write(&out_path, &e))?;
    let (mut rows, mut marked) = (0u64, 0u64);
    let mut last: Option<i64> = None;
    let mut n = 0u64;
    while let Some(raw) = t.raw()? {
        n += 1;
        let line = raw.strip_suffix('\n').unwrap_or(raw);
        let mut c: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if c.len() != TABLE_HEADER.len() {
            return die(format!(
                "{rows_path}: a row has {} columns, want {}: {}",
                c.len(),
                TABLE_HEADER.len(),
                head80(raw)
            ));
        }
        let pos = need_int(&c[1])?;
        // The rows come sorted (Picard sorts its output), so a position
        // seen before is the one just before: a repeat keeps its first row.
        match last {
            Some(l) if pos == l => continue,
            Some(l) if pos < l => {
                return usage(format!(
                    "{rows_path}: line {n}: position {pos} comes after {l}; the lifted rows must be sorted by position"
                ));
            }
            _ => {}
        }
        last = Some(pos);
        let p = pos - 1;
        let i = iv.partition_point(|&(s, _)| s <= p);
        let inside = i > 0 && p < iv[i - 1].1;
        if !inside {
            let mut hit = false;
            for cell in &mut c[4..7] {
                if let Some(v) = ["0", "1", "2"].iter().position(|x| x == cell) {
                    *cell = (v as u8 + MASKED).to_string();
                    hit = true;
                }
            }
            marked += u64::from(hit);
        }
        out.write_all((c.join("\t") + "\n").as_bytes())
            .map_err(|e| cannot_write(&out_path, &e))?;
        rows += 1;
    }
    out.commit().map_err(|e| cannot_write(&out_path, &e))?;
    let _ = writeln!(
        err,
        "{chrom}: {rows} rows, {marked} outside the strict mask and marked masked"
    );
    Ok(())
}

// --- genotypes (step 43) ----------------------------------------------------------

/// The stand-in columns for each ALT allele count 0..2n: as many 2s and
/// one 1 as the count needs, then 0s.
fn pseudo_columns(n: usize) -> Vec<String> {
    (0..=2 * n)
        .map(|ac| {
            let mut g = vec!["2"; ac / 2];
            if ac % 2 == 1 {
                g.push("1");
            }
            g.resize(n, "0");
            g.join("\t")
        })
        .collect()
}

/// A call of the sample's VCF: REF, ALT (`None` for an insertion, a
/// deletion or a multi-allelic record) and the code (0, 1, 2, or 9 when it
/// cannot be read).
#[derive(Debug, Clone, PartialEq, Eq)]
struct Sample {
    r#ref: String,
    alt: Option<String>,
    code: u8,
}

/// What [`pick_call`] keeps.
#[derive(Debug, PartialEq, Eq)]
enum Picked {
    Call(Sample),
    Clash,
}

/// The call to keep when a position comes twice, as in a VCF split with
/// bcftools norm -m-. An SNV wins over an insertion or a multi-allelic
/// record, in either order. Of one SNV twice, a usable record wins over
/// code 9, and two different usable codes clash. Of two different SNVs,
/// the one that carries its ALT wins over a 0/0 or unusable one; both 0/0
/// is the reference; both carrying their ALT, or neither readable, clash
/// (unknown, 9).
fn pick_call(prev: Option<Sample>, call: Sample) -> Picked {
    let Some(prev) = prev else {
        return Picked::Call(call);
    };
    if prev.alt.is_none() || call.alt.is_none() {
        return Picked::Call(if prev.alt.is_none() && call.alt.is_some() {
            call
        } else {
            prev
        });
    }
    if call.alt == prev.alt {
        if call.code == 9 || call.code == prev.code {
            return Picked::Call(prev);
        }
        return if prev.code == 9 {
            Picked::Call(call)
        } else {
            Picked::Clash
        };
    }
    let carries = |c: u8| c == 1 || c == 2;
    let (has_prev, has_call) = (carries(prev.code), carries(call.code));
    if has_prev != has_call {
        return Picked::Call(if has_prev { prev } else { call });
    }
    if !has_prev && prev.code == 0 && call.code == 0 {
        return Picked::Call(Sample {
            r#ref: prev.r#ref,
            alt: None,
            code: 0,
        });
    }
    Picked::Clash
}

/// One record of the sample's VCF that `genotypes` reads.
struct VcfRec {
    chrom: String,
    pos: i64,
    call: Sample,
}

/// The biallelic SNVs (and the other records, as code 9) of a one-sample
/// VCF, in file order. A record whose FILTER is not PASS or `.` gives
/// code 9 (not read), so a DeepVariant RefCall is neither your allele nor
/// the reference.
struct Vcf {
    text: Text,
}

impl Vcf {
    fn next(&mut self) -> Out<Option<VcfRec>> {
        loop {
            let Some(raw) = self.text.raw()? else {
                return Ok(None);
            };
            if raw.starts_with("##") {
                continue;
            }
            let line = raw.strip_suffix('\n').unwrap_or(raw);
            let c: Vec<&str> = line.split('\t').collect();
            if raw.starts_with('#') {
                if c.len() != 10 {
                    let n = c.len().saturating_sub(9);
                    return die(format!(
                        "{} holds {n} samples; step 43 reads a one-sample VCF \
                         (keep one with bcftools view -s <sample>).",
                        self.text.path.display()
                    ));
                }
                continue;
            }
            if c.len() < 10 || c[3].chars().count() != 1 || c[4] == "." {
                continue;
            }
            if c[4].contains(',') || c[4].chars().count() != 1 {
                return Ok(Some(VcfRec {
                    chrom: c[0].to_owned(),
                    pos: need_int(c[1])?,
                    call: Sample {
                        r#ref: c[3].to_owned(),
                        alt: None,
                        code: 9,
                    },
                }));
            }
            let fmt: Vec<&str> = c[8].split(':').collect();
            let Some(k) = fmt.iter().position(|f| *f == "GT") else {
                continue;
            };
            let gt = c[9].split(':').nth(k).unwrap_or(".");
            let mut code = gt_code(gt).unwrap_or(9);
            if c[6] != "PASS" && c[6] != "." {
                code = 9;
            }
            return Ok(Some(VcfRec {
                chrom: c[0].to_owned(),
                pos: need_int(c[1])?,
                call: Sample {
                    r#ref: c[3].to_owned(),
                    alt: Some(c[4].to_owned()),
                    code,
                },
            }));
        }
    }
}

/// The sample's calls on one chromosome, one per position after
/// [`pick_call`], sorted by position: a compact array, not a map.
fn settle(mut recs: Vec<(i64, Sample)>) -> Vec<(i64, Sample)> {
    // A stable sort keeps the records of one position in file order.
    recs.sort_by_key(|(p, _)| *p);
    let mut out: Vec<(i64, Sample)> = Vec::new();
    let mut clashed = false;
    for (pos, call) in recs {
        match out.last_mut() {
            Some((p, kept)) if *p == pos => {
                if clashed {
                    continue;
                }
                let r#ref = call.r#ref.clone();
                match pick_call(Some(kept.clone()), call) {
                    Picked::Call(c) => *kept = c,
                    Picked::Clash => {
                        *kept = Sample {
                            r#ref,
                            alt: None,
                            code: 9,
                        };
                        clashed = true;
                    }
                }
            }
            _ => {
                clashed = false;
                out.push((pos, call));
            }
        }
    }
    out
}

/// `chrN.tsv.gz` with N from 1 to 22: the chromosome's name.
fn table_chrom(name: &str) -> Option<&str> {
    let chrom = name.strip_suffix(".tsv.gz")?;
    let num = chrom.strip_prefix("chr")?;
    let ok = match num.len() {
        1 => matches!(num.as_bytes()[0], b'1'..=b'9'),
        2 => {
            let (a, b) = (num.as_bytes()[0], num.as_bytes()[1]);
            (a == b'1' && b.is_ascii_digit()) || (a == b'2' && (b'0'..=b'2').contains(&b))
        }
        _ => false,
    };
    ok.then_some(chrom)
}

/// A file written whole through a temporary name.
fn put(path: &Path, data: &[u8]) -> Result<()> {
    let mut f = AtomicFile::create(path).map_err(|e| cannot_write(path, &e))?;
    f.write_all(data).map_err(|e| cannot_write(path, &e))?;
    f.commit().map_err(|e| cannot_write(path, &e))
}

fn cmd_genotypes(a: &Args, err: &mut dyn Write) -> Out<()> {
    let pseudo_n = a
        .int("--pseudo", PSEUDO)
        .map_err(|_| Stop(2, String::new()))?;
    if pseudo_n < 1 {
        return usage("--pseudo must be 1 or more");
    }
    let pop = a.get("--pop").unwrap_or("ALL");
    let Some(af_col) = POPS.iter().position(|p| *p == pop).map(|k| 7 + k) else {
        return usage(format!(
            "--pop must be one of {}, got {}",
            POPS.join(", "),
            repr(pop)
        ));
    };
    let table_dir = a.req("--table-dir");
    let rd = std::fs::read_dir(table_dir).map_err(|e| cannot_read(Path::new(table_dir), &e))?;
    let mut names: Vec<String> = Vec::new();
    for e in rd {
        let e = e.map_err(|e| cannot_read(Path::new(table_dir), &e))?;
        if let Ok(n) = e.file_name().into_string() {
            names.push(n);
        }
    }
    names.sort();
    let tables: Vec<(String, String)> = names
        .iter()
        .filter_map(|n| table_chrom(n).map(|c| (c.to_owned(), join(table_dir, n))))
        .collect();
    if tables.is_empty() {
        return die(format!(
            "{table_dir} holds no chrN.tsv.gz: install the archaic data with scripts/setup.sh --archaic"
        ));
    }
    let out_dir = a.req("--out-dir");
    std::fs::create_dir_all(out_dir).map_err(|e| cannot_write(Path::new(out_dir), &e))?;
    let n = pseudo_n as usize;
    let pseudo = pseudo_columns(n);
    let two_n = (2 * n) as f64;
    let mut head = String::from("chrom\tpos\tref\talt\tAltai\tVindija\tDenisova\tSAMPLE");
    let mut samples = String::from("SAMPLE\n");
    for i in 1..=n {
        let id = format!("P{}", zero_pad(i as i128, 3));
        head += &format!("\t{id}");
        samples += &format!("{id}\n");
    }
    head.push('\n');
    put(&Path::new(out_dir).join("samples.txt"), samples.as_bytes())?;

    let vcf_path = a.req("--vcf");
    let mut vcf = Vcf {
        text: Text::open(vcf_path)?,
    };
    let mut spans: Vec<Span> = Vec::new();
    let mut done: Vec<String> = Vec::new();
    let mut pending = vcf.next()?;
    while let Some(first) = pending.take() {
        let chrom = first.chrom.clone();
        if done.contains(&chrom) {
            return die(format!(
                "{vcf_path} is not sorted: {chrom} comes back after another chromosome"
            ));
        }
        done.push(chrom.clone());
        let table = tables.iter().find(|(c, _)| *c == chrom).map(|(_, t)| t);
        // The sample's calls on this chromosome; on a chromosome the table
        // does not hold they are read past.
        let mut recs = Vec::new();
        if table.is_some() {
            recs.push((first.pos, first.call));
        }
        loop {
            match vcf.next()? {
                Some(r) if r.chrom == chrom => {
                    if table.is_some() {
                        recs.push((r.pos, r.call));
                    }
                }
                other => {
                    pending = other;
                    break;
                }
            }
        }
        let Some(table) = table else {
            continue;
        };
        let calls = settle(recs);
        if let Some(span) = write_chrom(
            &chrom, table, &calls, &head, &pseudo, af_col, two_n, out_dir,
        )? {
            spans.push(span);
        }
    }
    let mut text =
        String::from("chrom\tfirst\tlast\tsites\tsites_with_your_alt\tsites_other_allele\n");
    for (c, f, l, n, carried, dropped) in &spans {
        text += &format!("{c}\t{f}\t{l}\t{n}\t{carried}\t{dropped}\n");
    }
    put(&Path::new(out_dir).join("spans.tsv"), text.as_bytes())?;
    if spans.is_empty() {
        let names: Vec<&str> = tables.iter().map(|(c, _)| c.as_str()).collect();
        let mut sorted = names.clone();
        sorted.sort_unstable();
        return die(format!(
            "{vcf_path} has no records on a chromosome of the archaic table ({})",
            sorted.join(", ")
        ));
    }
    let sites: u64 = spans.iter().map(|s| s.3).sum();
    let _ = writeln!(err, "{} chromosomes, {sites} sites", spans.len());
    Ok(())
}

/// A chromosome's row of `spans.tsv`: chromosome, first and last
/// position, sites, sites with your ALT, sites with another allele.
type Span = (String, i64, i64, u64, u64, u64);

/// IBDmix's genotype file of one chromosome, and each archaic genome's
/// mask; the chromosome's span row, or `None` when no row was written.
#[allow(clippy::too_many_arguments)]
fn write_chrom(
    chrom: &str,
    table: &str,
    calls: &[(i64, Sample)],
    head: &str,
    pseudo: &[String],
    af_col: usize,
    two_n: f64,
    out_dir: &str,
) -> Out<Option<Span>> {
    let num = &chrom[3..];
    let path = Path::new(out_dir).join(format!("{chrom}.gt.gz"));
    let file = AtomicFile::create(&path).map_err(|e| cannot_write(&path, &e))?;
    let mut gz = GzEncoder::new(file, Compression::new(1));
    let werr = |e: io::Error| Stop::from(cannot_write(&path, &e));
    gz.write_all(head.as_bytes()).map_err(werr)?;
    let mut t = Text::open(table)?;
    let header = t.next()?.unwrap_or("");
    if header.split('\t').ne(TABLE_HEADER.iter().copied()) {
        return die(format!(
            "{table} does not start with the table header ({})",
            TABLE_HEADER.join(" ")
        ));
    }
    let (mut n, mut first, mut last, mut dropped, mut carried) = (0u64, 0i64, 0i64, 0u64, 0u64);
    // Per archaic genome, the runs of written rows whose call is masked:
    // [first, last] positions, closed at the next unmasked row.
    let mut runs: [Vec<(i64, i64)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut open_run: [Option<(i64, i64)>; 3] = [None; 3];
    let mut row = String::new();
    while let Some(line) = t.next()? {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 7 || c.len() <= af_col {
            return die(format!(
                "{table}: a row has {} columns, want {}",
                c.len(),
                TABLE_HEADER.len()
            ));
        }
        let pos = need_int(c[1])?;
        let (r, alt) = (c[2], c[3]);
        let mut arch = [0u8; 3];
        for k in 0..3 {
            let v = c[4 + k];
            arch[k] = match v {
                "0" | "1" | "2" | "9" | "10" | "11" | "12" => v.parse().expect("a code"),
                _ => {
                    return die(format!(
                        "{table}: {chrom}:{pos} has {} code {}; the table holds 0, 1, 2, 9 or {MASKED} to {}",
                        ARCHAICS[k],
                        repr(v),
                        MASKED + 2
                    ));
                }
            };
        }
        let code = match calls.binary_search_by_key(&pos, |(p, _)| *p) {
            Err(_) => 0,
            Ok(i) => {
                let call = &calls[i].1;
                if call.r#ref != r || call.alt.as_ref().is_some_and(|x| x != alt) {
                    // Another allele than the row's: IBDmix's generate_gt
                    // drops the line.
                    dropped += 1;
                    continue;
                }
                call.code
            }
        };
        if code == 1 || code == 2 {
            carried += 1;
        }
        let ac = match pycompat::num::parse_float(c[af_col]) {
            Err(_) => 0.0,
            Ok(x) => match pycompat::fmt::round_int(x * two_n) {
                Ok(v) => v,
                Err(pycompat::fmt::RoundError::Nan) => 0.0,
                Err(e) => return die(format!("{table}: {chrom}:{pos}: {e}")),
            },
        };
        let ac = ac.clamp(0.0, two_n) as usize;
        for k in 0..3 {
            if arch[k] >= MASKED {
                match &mut open_run[k] {
                    None => open_run[k] = Some((pos, pos)),
                    Some(run) => run.1 = pos,
                }
            } else if let Some(run) = open_run[k].take() {
                runs[k].push(run);
            }
        }
        row.clear();
        row += &format!("{num}\t{pos}\t{r}\t{alt}");
        for v in arch {
            let v = if v >= MASKED { v - MASKED } else { v };
            row += &format!("\t{v}");
        }
        row += &format!("\t{code}\t{}\n", pseudo[ac]);
        gz.write_all(row.as_bytes()).map_err(werr)?;
        if n == 0 {
            first = pos;
        }
        last = pos;
        n += 1;
    }
    if n == 0 {
        drop(gz);
        match std::fs::remove_file(&path) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => {
                return Err(Stop::from(cannot_write(&path, &e)));
            }
            _ => {}
        }
        return Ok(None);
    }
    gz.finish().and_then(AtomicFile::commit).map_err(werr)?;
    // IBDmix's mask (-r): the chromosome as in the genotype file, and a
    // line covers the positions start + 1 to end, as in a BED.
    for k in 0..3 {
        if let Some(run) = open_run[k] {
            runs[k].push(run);
        }
        let text: String = runs[k]
            .iter()
            .map(|(s, e)| format!("{num}\t{}\t{e}\n", s - 1))
            .collect();
        put(
            &Path::new(out_dir).join(format!("{chrom}.{}.mask.bed", ARCHAICS[k])),
            text.as_bytes(),
        )?;
    }
    Ok(Some((chrom.to_owned(), first, last, n, carried, dropped)))
}

// --- report (step 43) -------------------------------------------------------------

/// A segment: chromosome, start, end, LOD.
type Seg = (String, i64, i64, f64);

/// IBDmix's output rows for SAMPLE, kept at LOD 4 and 50 kb or longer.
fn read_segments(path: &str) -> Out<Vec<Seg>> {
    let mut t = Text::open(path)?;
    let head: Vec<String> = t
        .next()?
        .unwrap_or("")
        .split('\t')
        .map(str::to_owned)
        .collect();
    let first5: Vec<&str> = head.iter().take(5).map(String::as_str).collect();
    if first5 != ["ID", "chrom", "start", "end", "slod"] {
        return die(format!(
            "{path} is not IBDmix output (header {})",
            first5.join(" ")
        ));
    }
    let mut segs = Vec::new();
    while let Some(line) = t.next()? {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 || c[0] != "SAMPLE" {
            continue;
        }
        let (s, e, lod) = (need_int(c[2])?, need_int(c[3])?, need_float(c[4])?);
        if lod >= MIN_LOD && e - s >= MIN_LENGTH {
            let chrom = if c[1].starts_with("chr") {
                c[1].to_owned()
            } else {
                format!("chr{}", c[1])
            };
            segs.push((chrom, s, e, lod));
        }
    }
    Ok(segs)
}

/// The order of chromosomes in the report: the number, or 99 for any
/// other name (so chrX and chrY sort together, by position).
fn chrom_key(c: &str) -> i128 {
    let c = c.strip_prefix("chr").unwrap_or(c);
    if pycompat::text::isdigit(c) {
        pycompat::num::parse_int(c).unwrap_or(99)
    } else {
        99
    }
}

/// A merged segment: chromosome, start, end, the highest LOD, and which
/// archaic genomes it holds.
type Merged = (String, i64, i64, f64, Vec<&'static str>);

/// Overlapping or touching segments of one class as one.
fn merge(mut segs: Vec<(Seg, &'static str)>) -> Vec<Merged> {
    segs.sort_by(|(a, _), (b, _)| (chrom_key(&a.0), a.1, a.2).cmp(&(chrom_key(&b.0), b.1, b.2)));
    let mut out: Vec<Merged> = Vec::new();
    for ((chrom, start, end, lod), who) in segs {
        match out.last_mut() {
            Some(m) if m.0 == chrom && start <= m.2 => {
                m.2 = m.2.max(end);
                m.3 = py_max(m.3, lod);
                if !m.4.contains(&who) {
                    m.4.push(who);
                }
            }
            _ => out.push((chrom, start, end, lod, vec![who])),
        }
    }
    out
}

fn overlaps(chrom: &str, s: i64, e: i64, others: &[Merged]) -> bool {
    others.iter().any(|o| o.0 == chrom && o.1 < e && s < o.2)
}

/// A row of the report: chromosome, start, end, ancestry, matched, LOD.
type Row = (String, i64, i64, &'static str, Vec<&'static str>, f64);

/// Neanderthal: the union of the Altai and Vindija segments. Denisovan:
/// the Denisova segments that overlap no Neanderthal one (a stretch both
/// match is counted as Neanderthal, the more common source outside
/// Oceania).
fn classify(per: &[Vec<Seg>; 3]) -> Vec<Row> {
    let nea = merge(
        NEANDERTHALS
            .iter()
            .enumerate()
            .flat_map(|(k, n)| per[k].iter().map(move |s| (s.clone(), *n)))
            .collect(),
    );
    let den = merge(per[2].iter().map(|s| (s.clone(), "Denisova")).collect());
    let mut rows: Vec<Row> = Vec::new();
    for (c, s, e, l, mut who) in nea.iter().cloned() {
        who.sort_unstable();
        if overlaps(&c, s, e, &den) {
            who.push("Denisova");
        }
        rows.push((c, s, e, "Neanderthal", who, l));
    }
    for (c, s, e, l, _) in den.iter().cloned() {
        if !overlaps(&c, s, e, &nea) {
            rows.push((c, s, e, "Denisovan", vec!["Denisova"], l));
        }
    }
    rows.sort_by(|a, b| (chrom_key(&a.0), a.1).cmp(&(chrom_key(&b.0), b.1)));
    rows
}

/// Genes per chromosome: (start, end, name), sorted.
type Genes = std::collections::HashMap<String, Vec<(i64, i64, String)>>;

/// Protein-coding genes of a GENCODE gene GTF, per chromosome, sorted by
/// (start, end, name).
fn read_genes(path: &str) -> Out<Genes> {
    let mut genes = Genes::new();
    let mut t = Text::open(path)?;
    while let Some(line) = t.next()? {
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 9 || c[2] != "gene" || !c[8].contains("gene_type \"protein_coding\"") {
            continue;
        }
        if let Some(name) = gene_name(c[8]) {
            let (s, e) = (need_int(c[3])?, need_int(c[4])?);
            genes
                .entry(c[0].to_owned())
                .or_default()
                .push((s, e, name.to_owned()));
        }
    }
    for v in genes.values_mut() {
        v.sort();
    }
    Ok(genes)
}

/// The first `gene_name "X"` with a non-empty X (the regex
/// `gene_name "([^"]+)"`).
fn gene_name(attrs: &str) -> Option<&str> {
    const KEY: &str = "gene_name \"";
    attrs.match_indices(KEY).find_map(|(i, _)| {
        let rest = &attrs[i + KEY.len()..];
        let end = rest.find('"')?;
        (end > 0).then(|| &rest[..end])
    })
}

/// The names of the genes overlapping `[start, end)`, once each, in order.
/// IBDmix's start is the segment's first site and its end the next site,
/// past the segment; a GTF gene is `[s, e]`, 1-based.
fn genes_in<'a>(genes: &'a Genes, chrom: &str, start: i64, end: i64) -> Vec<&'a str> {
    let name = if chrom.starts_with("chr") {
        chrom.to_owned()
    } else {
        format!("chr{chrom}")
    };
    let mut seen: Vec<&str> = Vec::new();
    for (s, e, n) in genes.get(&name).map_or(&[][..], Vec::as_slice) {
        if *s >= end {
            break; // sorted by start: no later gene starts before the end
        }
        if *e >= start && !seen.contains(&n.as_str()) {
            seen.push(n);
        }
    }
    seen
}

fn cmd_report(a: &Args, out: &mut dyn Write) -> Out<()> {
    let spans_path = a.req("--spans");
    let mut t = Text::open(spans_path)?;
    let head = t.next()?.unwrap_or("").to_owned();
    if head
        .split('\t')
        .take(4)
        .ne(["chrom", "first", "last", "sites"])
    {
        return die(format!(
            "{spans_path} is not the spans table of the genotypes step"
        ));
    }
    let mut spans: Vec<(String, i64, i64, i64)> = Vec::new();
    while let Some(line) = t.next()? {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 4 {
            return die(format!(
                "{spans_path}: a row has {} columns, want 4",
                c.len()
            ));
        }
        spans.push((
            c[0].to_owned(),
            need_int(c[1])?,
            need_int(c[2])?,
            need_int(c[3])?,
        ));
    }
    if spans.is_empty() {
        return die(format!("{spans_path} lists no chromosome"));
    }
    let spans_bp: i64 = spans.iter().map(|(_, f, l, _)| l - f + 1).sum();
    let ibdmix = a.req("--ibdmix-dir");
    let mut per: [Vec<Seg>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for (chrom, ..) in &spans {
        for (k, n) in ARCHAICS.iter().enumerate() {
            let p = join(ibdmix, &format!("{chrom}.{n}.txt"));
            if !Path::new(&p).is_file() {
                return die(format!("no IBDmix output {p}"));
            }
            per[k].extend(read_segments(&p)?);
        }
    }
    let rows = classify(&per);
    let genes_arg = a.get("--genes").unwrap_or("");
    let genes = if genes_arg.is_empty() {
        Default::default()
    } else {
        read_genes(genes_arg)?
    };
    let genes_of = |c: &str, s: i64, e: i64| -> Vec<&str> {
        if genes_arg.is_empty() {
            Vec::new()
        } else {
            genes_in(&genes, c, s, e)
        }
    };
    let sample = a.req("--sample");
    let pop = a.get("--pop").unwrap_or("ALL");
    let out_dir = a.req("--out-dir");
    std::fs::create_dir_all(out_dir).map_err(|e| cannot_write(Path::new(out_dir), &e))?;

    let seg_name = format!("{sample}_archaic_segments.tsv");
    let mut text = format!("{SEGMENT_HEADER}\n");
    for (c, s, e, anc, who, lod) in &rows {
        let g = genes_of(c, *s, *e);
        text += &format!(
            "{c}\t{s}\t{e}\t{}\t{anc}\t{}\t{}\t{}\n",
            e - s,
            who.join(","),
            fixed(*lod, 1),
            if g.is_empty() {
                ".".to_owned()
            } else {
                g.join(",")
            }
        );
    }
    put(&Path::new(out_dir).join(&seg_name), text.as_bytes())?;

    // The summary of each ancestry.
    let mut summary: Vec<(String, String)> = Vec::new();
    for anc in ["Neanderthal", "Denisovan"] {
        let r: Vec<&Row> = rows.iter().filter(|x| x.3 == anc).collect();
        let bp: i64 = r.iter().map(|x| x.2 - x.1).sum();
        let key = anc.to_lowercase();
        summary.push((format!("{key}_segments"), r.len().to_string()));
        summary.push((format!("{key}_bp"), bp.to_string()));
        summary.push((
            format!("{key}_percent"),
            if spans_bp != 0 {
                fixed(100.0 * bp as f64 / spans_bp as f64, 2)
            } else {
                "0.00".to_owned()
            },
        ));
        summary.push((
            format!("{key}_longest_bp"),
            r.iter().map(|x| x.2 - x.1).max().unwrap_or(0).to_string(),
        ));
    }
    let sv = |k: &str| {
        summary
            .iter()
            .find(|(key, _)| key == k)
            .map(|(_, v)| v.as_str())
            .expect("a summary key")
    };
    let mut longest: Vec<&Row> = rows.iter().collect();
    longest.sort_by(|a, b| {
        (-(a.2 - a.1), chrom_key(&a.0), a.1).cmp(&(-(b.2 - b.1), chrom_key(&b.0), b.1))
    });
    longest.truncate(TOP);
    let sites: i64 = spans.iter().map(|s| s.3).sum();
    let mut text = String::from("key\tvalue\n");
    let meta: [(&str, String); 10] = [
        ("sample", sample.to_owned()),
        ("method", "IBDmix".to_owned()),
        (
            "archaic_genomes",
            "Altai Neanderthal, Vindija 33.19 Neanderthal, Altai Denisovan".to_owned(),
        ),
        ("population", pop.to_owned()),
        ("min_lod", general(MIN_LOD, 6)),
        ("min_length_bp", MIN_LENGTH.to_string()),
        ("analysed_bp", spans_bp.to_string()),
        ("analysed_chromosomes", spans.len().to_string()),
        ("analysed_sites", sites.to_string()),
        (
            "genes",
            if genes_arg.is_empty() {
                "none (no gene annotation installed)"
            } else {
                "GENCODE protein-coding"
            }
            .to_owned(),
        ),
    ];
    for (k, v) in meta.iter() {
        text += &format!("{k}\t{v}\n");
    }
    for (k, v) in &summary {
        text += &format!("{k}\t{v}\n");
    }
    let genes_text = |g: &[&str]| {
        if g.is_empty() {
            "none".to_owned()
        } else {
            g.join(", ")
        }
    };
    for (i, (c, s, e, anc, who, _)) in longest.iter().enumerate() {
        let g = genes_of(c, *s, *e);
        text += &format!(
            "longest_{}\t{c}:{s}-{e} {anc} {} Mb ({}); genes: {}\n",
            i + 1,
            fixed((e - s) as f64 / 1e6, 2),
            who.join("/"),
            genes_text(&g)
        );
    }
    put(
        &Path::new(out_dir).join(format!("{sample}_archaic.tsv")),
        text.as_bytes(),
    )?;
    let mb = |k: &str| fixed(need_int(sv(k)).unwrap_or(0) as f64 / 1e6, 1);
    let mut lines = vec![
        format!(
            "Archaic ancestry of {sample} (IBDmix against the Altai and Vindija Neanderthals and the Altai Denisovan)"
        ),
        format!(
            "  Analysed: {} autosomes, {} Mb, {sites} sites; modern frequencies: 1000 Genomes {pop}",
            spans.len(),
            fixed(spans_bp as f64 / 1e6, 1)
        ),
        format!(
            "  Neanderthal: {}% ({} Mb in {} segments)",
            sv("neanderthal_percent"),
            mb("neanderthal_bp"),
            sv("neanderthal_segments")
        ),
        format!(
            "  Denisovan:   {}% ({} Mb in {} segments)",
            sv("denisovan_percent"),
            mb("denisovan_bp"),
            sv("denisovan_segments")
        ),
    ];
    if !longest.is_empty() {
        lines.push(format!("  The {} longest segments:", longest.len()));
        for (c, s, e, anc, who, lod) in &longest {
            let g = genes_of(c, *s, *e);
            lines.push(format!(
                "    {c}:{s}-{e}  {} Mb  {anc} ({}, LOD {})  genes: {}",
                fixed((e - s) as f64 / 1e6, 2),
                who.join("/"),
                fixed(*lod, 0),
                genes_text(&g)
            ));
        }
    }
    lines.push(format!("  {CAVEAT}"));
    lines.push(format!("  Every segment: {seg_name}"));
    let body = lines.join("\n") + "\n";
    put(
        &Path::new(out_dir).join(format!("{sample}_archaic.txt")),
        body.as_bytes(),
    )?;
    let _ = out.write_all(body.as_bytes());
    Ok(())
}

// --- the command line -------------------------------------------------------------

const USAGE: &str = "usage: og archaic {table,strict,genotypes,report} ...\n";

const TABLE: Spec = Spec {
    verb: "table",
    usage: "usage: og archaic table --chrom CHROM --modern MODERN [--archaic NAME=FILE] [--mask NAME=BED] [--exclude BED] --out OUT\n",
    options: &[
        "--chrom",
        "--modern",
        "--archaic",
        "--mask",
        "--exclude",
        "--out",
    ],
    required: &["--chrom", "--modern", "--out"],
};
const STRICT: Spec = Spec {
    verb: "strict",
    usage: "usage: og archaic strict --chrom CHROM --strict STRICT --rows ROWS --out OUT\n",
    options: &["--chrom", "--strict", "--rows", "--out"],
    required: &["--chrom", "--strict", "--rows", "--out"],
};
const GENOTYPES: Spec = Spec {
    verb: "genotypes",
    usage: "usage: og archaic genotypes --table-dir TABLE_DIR --vcf VCF [--pop POP] [--pseudo PSEUDO] --out-dir OUT_DIR\n",
    options: &["--table-dir", "--vcf", "--pop", "--pseudo", "--out-dir"],
    required: &["--table-dir", "--vcf", "--out-dir"],
};
const REPORT: Spec = Spec {
    verb: "report",
    usage: "usage: og archaic report --sample SAMPLE --ibdmix-dir IBDMIX_DIR --spans SPANS [--genes GENES] [--pop POP] --out-dir OUT_DIR\n",
    options: &[
        "--sample",
        "--ibdmix-dir",
        "--spans",
        "--genes",
        "--pop",
        "--out-dir",
    ],
    required: &["--sample", "--ibdmix-dir", "--spans", "--out-dir"],
};

/// `og archaic <verb> [args]`: the exit code. Results go to `out`,
/// messages to `err`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let Some(verb) = args.first() else {
        let _ = write!(err, "{USAGE}");
        let _ = writeln!(
            err,
            "og archaic: error: the following arguments are required: verb"
        );
        return 2;
    };
    let spec = match verb.as_str() {
        "table" => &TABLE,
        "strict" => &STRICT,
        "genotypes" => &GENOTYPES,
        "report" => &REPORT,
        "-h" | "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        other => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(
                err,
                "og archaic: error: argument verb: invalid choice: {}",
                repr(other)
            );
            return 2;
        }
    };
    let a = match Args::parse(spec, &args[1..]).and_then(|a| {
        if spec.verb == "genotypes" {
            a.int("--pseudo", PSEUDO)?;
        }
        Ok(a)
    }) {
        Ok(a) => a,
        Err(ArgError::Help) => {
            let _ = write!(out, "{}", spec.usage);
            return 0;
        }
        Err(ArgError::Usage(msg)) => {
            let _ = write!(err, "{}", spec.usage);
            let _ = writeln!(err, "og archaic {}: error: {msg}", spec.verb);
            return 2;
        }
    };
    let done = match spec.verb {
        "table" => cmd_table(&a, err),
        "strict" => cmd_strict(&a, err),
        "genotypes" => cmd_genotypes(&a, err),
        _ => cmd_report(&a, out),
    };
    match done {
        Ok(()) => 0,
        Err(Stop(code, msg)) => {
            let _ = out.flush();
            let _ = writeln!(err, "ERROR: {msg}");
            code
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(r: &str, alt: Option<&str>, code: u8) -> Sample {
        Sample {
            r#ref: r.into(),
            alt: alt.map(Into::into),
            code,
        }
    }

    fn code(p: Picked) -> u8 {
        match p {
            Picked::Call(c) => c.code,
            Picked::Clash => 9,
        }
    }

    #[test]
    fn a_mask_holds_its_bed_bases_only() {
        // BED 99-100 and 200-300 are the 1-based bases 100 and 201 to 300.
        let mut m = Mask {
            iv: vec![(99, 100), (200, 300)],
            i: 0,
        };
        let got: Vec<i64> = [99, 100, 101, 200, 201, 300, 301]
            .into_iter()
            .filter(|&p| m.holds(p))
            .collect();
        assert_eq!(got, [100, 201, 300]);
    }

    #[test]
    fn gt_codes() {
        assert_eq!(gt_code("0/1"), Some(1));
        assert_eq!(gt_code("1|1"), Some(2));
        assert_eq!(gt_code("0/0"), Some(0));
        assert_eq!(gt_code("./."), None);
        assert_eq!(gt_code("1/2"), None);
        assert_eq!(gt_code("1"), None);
        assert_eq!(gt_code("0/1/1"), None);
    }

    #[test]
    fn pick_call_is_symmetric_in_its_code() {
        let pairs = [
            (s("A", Some("G"), 1), s("A", None, 9)),
            (s("A", Some("G"), 1), s("A", Some("T"), 0)),
            (s("A", Some("G"), 9), s("A", Some("G"), 2)),
            (s("A", Some("G"), 0), s("A", Some("T"), 0)),
            (s("A", Some("G"), 1), s("A", Some("T"), 2)),
            (s("A", Some("G"), 9), s("A", Some("T"), 0)),
            (s("A", Some("G"), 1), s("A", Some("G"), 2)),
        ];
        for (x, y) in pairs {
            assert_eq!(
                code(pick_call(Some(x.clone()), y.clone())),
                code(pick_call(Some(y), x))
            );
        }
        assert_eq!(
            pick_call(Some(s("A", Some("G"), 0)), s("A", Some("T"), 0)),
            Picked::Call(s("A", None, 0))
        );
    }

    #[test]
    fn a_clash_ends_its_position() {
        let recs = vec![
            (4000, s("T", Some("C"), 1)),
            (4000, s("T", Some("G"), 1)),
            (4000, s("T", Some("C"), 1)),
            (10, s("A", Some("G"), 2)),
        ];
        assert_eq!(
            settle(recs),
            vec![(10, s("A", Some("G"), 2)), (4000, s("T", None, 9))]
        );
    }

    #[test]
    fn pseudo_columns_carry_the_count() {
        let p = pseudo_columns(2);
        assert_eq!(p, ["0\t0", "1\t0", "2\t0", "2\t1", "2\t2"]);
    }

    #[test]
    fn table_names() {
        assert_eq!(table_chrom("chr1.tsv.gz"), Some("chr1"));
        assert_eq!(table_chrom("chr22.tsv.gz"), Some("chr22"));
        assert_eq!(table_chrom("chr19.tsv.gz"), Some("chr19"));
        assert_eq!(table_chrom("chr23.tsv.gz"), None);
        assert_eq!(table_chrom("chr0.tsv.gz"), None);
        assert_eq!(table_chrom("chrX.tsv.gz"), None);
        assert_eq!(table_chrom("chr01.tsv.gz"), None);
        assert_eq!(table_chrom("chr1.tsv"), None);
    }

    #[test]
    fn frequencies_clamp_like_python() {
        let afs: Vec<String> = ["1e-05", "nan", "inf", "-0.0", "x", "1.5"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let (_, _, _, clean) = table_row(
            Some(("A", "G", &afs)),
            [Some(("A", "G", "0/1")), None, None],
        )
        .unwrap();
        assert_eq!(clean, ["1e-05", "nan", "1", "-0", "0", "1"]);
    }

    #[test]
    fn genes_touching_the_ends() {
        let mut genes = std::collections::HashMap::new();
        let mut v = vec![
            (90, 100, "ENDS_AT_START".to_owned()),
            (100, 120, "INSIDE".to_owned()),
            (200, 210, "STARTS_AT_END".to_owned()),
            (50, 99, "BEFORE".to_owned()),
        ];
        v.sort();
        genes.insert("chr22".to_owned(), v);
        assert_eq!(
            genes_in(&genes, "22", 100, 200),
            ["ENDS_AT_START", "INSIDE"]
        );
    }

    #[test]
    fn gene_names() {
        assert_eq!(gene_name("gene_id \"G\"; gene_name \"ABC\";"), Some("ABC"));
        assert_eq!(gene_name("gene_name \"\"; gene_name \"B\";"), Some("B"));
        assert_eq!(gene_name("gene_name \"open"), None);
    }

    #[test]
    fn merge_keeps_the_highest_lod() {
        let m = merge(vec![
            (("chr22".into(), 1_000_000, 1_300_000, 40.0), "Altai"),
            (("chr22".into(), 1_200_000, 1_500_000, 50.0), "Vindija"),
            (("chr22".into(), 1_500_000, 1_600_000, 7.0), "Altai"),
        ]);
        assert_eq!(m.len(), 1);
        assert_eq!((m[0].1, m[0].2, m[0].3), (1_000_000, 1_600_000, 50.0));
        assert_eq!(m[0].4, ["Altai", "Vindija"]);
    }
}
