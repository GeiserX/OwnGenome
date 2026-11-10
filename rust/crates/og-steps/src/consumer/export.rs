//! `og consumer export`: every site of the table genotyped from the gVCF,
//! written in the 23andMe layout and the AncestryDNA layout, each zipped,
//! with a summary table.
//!
//! A site inside a reference block (GT 0/0) is homozygous reference; a site
//! with no gVCF record, a no-call, a filtered call or one inside a deletion
//! is `--` (AncestryDNA: `0 0`). Both files are read back with the strict
//! reader ([`super::check`]) before they are zipped. Exit 1 when a file
//! fails it, or when the gVCF's reference bases disagree with the site table
//! (the gVCF is not on GRCh38, or the table is wrong).
//!
//! The gVCF is read once, BGZF on many cores; only the lines near a site are
//! parsed. A site takes 24 bytes, so the whole chip list fits in about 15 MB.
//!
//! Deliberately unlike `bin/export_consumer.py`:
//! - the zips are the same bytes for the same input: each member carries
//!   the date of the header, not the time of the run, and the date is
//!   `--date`, else `SOURCE_DATE_EPOCH`, else today (UTC);
//! - the strict reader runs on the text before it is written, so no `.txt`
//!   file is left beside the zips (it still is when the reader fails);
//! - a call of more than two alleles stops the export, naming the site,
//!   where the Python wrote a file that then failed its own strict reader;
//! - a site row whose positions are not whole numbers from 0 to 2^32 - 1,
//!   or whose GRCh38 base is not one character, is refused as not a site
//!   row.

use std::fmt;
use std::io::{BufRead, Write};
use std::num::NonZero;
use std::time::{SystemTime, UNIX_EPOCH};

use og_io::{AtomicFile, Lines, Utf8};
use pycompat::num::parse_int;
use rustc_hash::FxHashSet;

use super::check::{self, Layout};
use super::zip::{self, DosDate};
use super::{
    Args, CHROMS, Fail, MT, Usage, X, Y, basename, chrom_index, gvcf_contig, head, path_join,
    py_repr, py_repr_list, report,
};

pub const DOC_URL: &str =
    "https://github.com/GeiserX/OwnGenome/blob/main/docs/40-export-consumer.md";
const SITES_HEADER: [&str; 7] = [
    "rsid", "chrom", "pos37", "ref37", "pos38", "ref38", "strand",
];
/// The X pseudoautosomal regions on GRCh37, 1-based inclusive (GRC).
const PAR37: [(u32, u32); 2] = [(60001, 2_699_520), (154_931_044, 155_260_560)];
/// ALT alleles that make a record a reference block.
const SYMBOLIC: [&[u8]; 3] = [b"<*>", b"<NON_REF>", b"."];

/// An rsID as the table spells it: `rs<n>` and `i<n>` without leading
/// zeros are kept as numbers; anything else as text in [`Table::names`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RsId {
    Rs(u32),
    I(u32),
    Name(u32),
}

impl RsId {
    fn new(s: &str, names: &mut Vec<Box<str>>) -> RsId {
        let number = |d: &str| {
            (!d.is_empty()
                && d.bytes().all(|c| c.is_ascii_digit())
                && (d == "0" || !d.starts_with('0')))
            .then(|| d.parse::<u32>().ok())
            .flatten()
        };
        if let Some(n) = s.strip_prefix("rs").and_then(number) {
            return RsId::Rs(n);
        }
        if let Some(n) = s.strip_prefix('i').and_then(number) {
            return RsId::I(n);
        }
        names.push(s.into());
        RsId::Name(u32::try_from(names.len() - 1).expect("fewer than 2^32 names"))
    }
}

/// Why a site has no genotype; [`fmt::Display`] is the reason the summary
/// counts it under.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Why {
    NoGenotype,
    NoCallInGvcf,
    Unreadable,
    InsideDeletion,
    RefMismatch,
    Filtered,
    Symbolic,
    NoRecord,
    ChrYFemale,
    HaploidHet,
}

impl fmt::Display for Why {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Why::NoGenotype => "no genotype",
            Why::NoCallInGvcf => "no-call in the gVCF",
            Why::Unreadable => "unreadable genotype",
            Why::InsideDeletion => "inside a deletion or a longer allele",
            Why::RefMismatch => "reference base differs from the site table",
            Why::Filtered => "filtered call",
            Why::Symbolic => "symbolic allele",
            Why::NoRecord => "no gVCF record",
            Why::ChrYFemale => "chrY site of a female sample",
            Why::HaploidHet => {
                "two alleles on a haploid site (chrM heteroplasmy, or male chrX/chrY)"
            }
        })
    }
}

/// The genotype of a site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Call {
    /// Not reached yet.
    Unset,
    No(Why),
    /// One or two bases, GRCh38 + strand (after [`finalise`]: the letters
    /// written, GRCh37 + strand, sorted).
    Bases {
        n: u8,
        b: [u8; 2],
    },
    /// A genotype of more than two alleles (never in a DeepVariant gVCF).
    Many(u8),
}

/// One site of the table.
#[derive(Debug, Clone, Copy)]
struct Site {
    rsid: RsId,
    pos37: u32,
    pos38: u32,
    chrom: u8,
    ref38: u8,
    minus: bool,
    call: Call,
}

/// The site table: its `## key: value` lines and its sites in file order.
struct Table {
    meta: Vec<(String, String)>,
    sites: Vec<Site>,
    names: Vec<Box<str>>,
}

impl Table {
    fn meta(&self, key: &str) -> Option<&str> {
        self.meta
            .iter()
            .rev()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    fn rsid(&self, id: RsId, buf: &mut String) {
        buf.clear();
        match id {
            RsId::Rs(n) => {
                buf.push_str("rs");
                buf.push_str(&n.to_string());
            }
            RsId::I(n) => {
                buf.push('i');
                buf.push_str(&n.to_string());
            }
            RsId::Name(i) => buf.push_str(&self.names[i as usize]),
        }
    }

    /// The order of two rsIDs as strings.
    fn rsid_cmp(&self, a: RsId, b: RsId) -> std::cmp::Ordering {
        if a == b {
            return std::cmp::Ordering::Equal;
        }
        let (mut x, mut y) = (String::new(), String::new());
        self.rsid(a, &mut x);
        self.rsid(b, &mut y);
        x.cmp(&y)
    }
}

fn read_table(path: &str) -> Result<Table, Fail> {
    let input = og_io::open(path).map_err(|e| Fail::io(path, &e))?;
    let mut lines = Lines::new(input, Utf8::Replace);
    let mut t = Table {
        meta: Vec::new(),
        sites: Vec::new(),
        names: Vec::new(),
    };
    let mut have_cols = false;
    let mut raw = String::new();
    let mut n = 0usize;
    while lines.read_line(&mut raw).map_err(|e| Fail::io(path, &e))? {
        n += 1;
        let line = raw.trim_end_matches(['\r', '\n']);
        if let Some(kv) = line.strip_prefix("## ") {
            let (k, v) = kv.split_once(": ").unwrap_or((kv, ""));
            t.meta.push((k.to_owned(), v.to_owned()));
            continue;
        }
        if let Some(cols) = line.strip_prefix('#') {
            let cols: Vec<&str> = cols.split('\t').collect();
            if cols != SITES_HEADER {
                return Err(Fail::input(format!(
                    "{path}: columns {}, want {}",
                    py_repr_list(&cols),
                    py_repr_list(&SITES_HEADER)
                )));
            }
            have_cols = true;
            continue;
        }
        if line.is_empty() {
            continue;
        }
        if !have_cols {
            return Err(Fail::input(format!(
                "{path}: no '#rsid ...' column line before the sites"
            )));
        }
        let bad = || {
            Fail::input(format!(
                "{path} line {n}: not a site row: {}",
                py_repr(head(line, 80))
            ))
        };
        let c: Vec<&str> = line.split('\t').collect();
        // `c[6] not in "+-"` is a substring test: "", "+", "-" and "+-" pass.
        if c.len() != 7 || !matches!(c[6], "" | "+" | "-" | "+-") {
            return Err(bad());
        }
        let chrom = chrom_index(c[1]).ok_or_else(bad)?;
        let pos = |s: &str| {
            parse_int(s)
                .ok()
                .and_then(|v| u32::try_from(v).ok())
                .ok_or_else(bad)
        };
        let (pos37, pos38) = (pos(c[2])?, pos(c[4])?);
        let [ref38] = c[5].as_bytes() else {
            return Err(bad());
        };
        let rsid = RsId::new(c[0], &mut t.names);
        t.sites.push(Site {
            rsid,
            pos37,
            pos38,
            chrom,
            ref38: *ref38,
            minus: c[6] == "-",
            call: Call::Unset,
        });
    }
    if t.sites.is_empty() {
        return Err(Fail::input(format!("{path} holds no site")));
    }
    t.sites.shrink_to_fit();
    Ok(t)
}

/// The pseudoautosomal regions of GRCh38 by contig, 1-based inclusive.
type Par = Vec<(String, i128, i128)>;

fn read_par(path: Option<&str>) -> Result<Par, Fail> {
    // `if not path`: no option and an empty one both mean no PAR.
    let Some(path) = path.filter(|p| !p.is_empty()) else {
        return Ok(Vec::new());
    };
    let text = std::fs::read(path).map_err(|e| Fail::io(path, &e))?;
    let text = String::from_utf8_lossy(&text);
    let mut par = Vec::new();
    for line in text.lines() {
        let c: Vec<&str> = pycompat::text::split(line, None, None);
        if c.len() >= 3 {
            let num = |s: &str| {
                parse_int(s)
                    .map_err(|_| Fail::data(format!("{path}: {} is not a number", py_repr(s))))
            };
            par.push((c[0].to_owned(), num(c[1])? + 1, num(c[2])?));
        }
    }
    Ok(par)
}

fn in_par(par: &Par, contig: &str, pos: u32) -> bool {
    let pos = i128::from(pos);
    par.iter()
        .any(|(c, s, e)| c == contig && *s <= pos && pos <= *e)
}

/// The counts the summary reports.
#[derive(Debug, Default)]
struct Tally {
    reference_checked: u64,
    reference_mismatch: u64,
    from_blocks: u64,
    from_variants: u64,
}

/// Python's `str.find(sub, start)` for one byte: -1 when absent.
fn find(s: &[u8], c: u8, start: isize) -> isize {
    let from = usize::try_from(start).unwrap_or(0);
    if from > s.len() {
        return -1;
    }
    memchr::memchr(c, &s[from..]).map_or(-1, |i| (from + i) as isize)
}

/// Python's `s[a:b]` with negative indices counted from the end.
fn slice(s: &[u8], a: isize, b: isize) -> &[u8] {
    let len = s.len() as isize;
    let fix = |i: isize| {
        if i < 0 { (len + i).max(0) } else { i.min(len) }
    };
    let (a, b) = (fix(a), fix(b));
    if a >= b {
        &[]
    } else {
        &s[a as usize..b as usize]
    }
}

/// Python's `int()` of a byte string, as i64.
fn int(s: &[u8]) -> Option<i64> {
    let v = parse_int(std::str::from_utf8(s).ok()?).ok()?;
    i64::try_from(v).ok()
}

/// One parsed gVCF record.
struct Rec<'a> {
    pos: i64,
    end: i64,
    ref_: &'a [u8],
    alts: Vec<&'a [u8]>,
    filter: &'a [u8],
    /// FORMAT and the first sample, when the line has both.
    sample: Option<(&'a [u8], &'a [u8])>,
}

impl Rec<'_> {
    fn is_block(&self) -> bool {
        self.alts
            .iter()
            .all(|a| SYMBOLIC.iter().any(|s| a.eq_ignore_ascii_case(s)))
    }

    /// The GT value: `.` when the record has none.
    fn gt(&self) -> Option<&[u8]> {
        let (fmt, smp) = self.sample?;
        let i = fmt.split(|&c| c == b':').position(|k| k == b"GT")?;
        smp.split(|&c| c == b':').nth(i)
    }
}

/// `re.split(r"[/|]", gt)`.
fn gt_parts(gt: &[u8]) -> impl Iterator<Item = &[u8]> {
    gt.split(|&c| c == b'/' || c == b'|')
}

/// The END of `(?:^|;)END=(\d+)` in INFO, as Python's `re.search` finds it.
fn info_end(info: &[u8]) -> Option<i64> {
    let mut at = 0;
    while let Some(i) = memchr::memmem::find(&info[at..], b"END=") {
        let k = at + i;
        if k == 0 || info[k - 1] == b';' {
            let digits = info[k + 4..]
                .iter()
                .take_while(|c| c.is_ascii_digit())
                .count();
            if digits > 0 {
                return int(&info[k + 4..k + 4 + digits]);
            }
        }
        at = k + 1;
    }
    None
}

/// `parse_record`: the record of a line (without its newline), or None for
/// fewer than eight columns.
fn parse_record(line: &[u8]) -> Result<Option<Rec<'_>>, String> {
    let c: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
    if c.len() < 8 {
        return Ok(None);
    }
    let pos = int(c[1]).ok_or_else(|| {
        format!(
            "position {} is not a number",
            py_repr(&String::from_utf8_lossy(c[1]))
        )
    })?;
    let mut end = pos + c[3].len() as i64 - 1;
    if let Some(e) = info_end(c[7]) {
        end = end.max(e);
    }
    Ok(Some(Rec {
        pos,
        end,
        ref_: c[3],
        alts: c[4].split(|&c| c == b',').collect(),
        filter: c[6],
        sample: (c.len() > 9).then(|| (c[8], c[9])),
    }))
}

/// `nonref_gt`: the record's genotype carries an allele other than REF.
fn nonref_gt(r: &Rec) -> bool {
    if r.sample.is_none() {
        return false;
    }
    r.gt()
        .is_some_and(|gt| gt_parts(gt).any(|i| !matches!(i, b"0" | b"." | b"")))
}

/// `call_from_record`: the bases at `pos` from one record, or why none.
fn call_from_record(r: &Rec, pos: i64, ref38: u8, tally: &mut Tally) -> Result<Call, String> {
    if r.sample.is_none() {
        return Ok(Call::No(Why::NoGenotype));
    }
    let gt = r.gt().unwrap_or(b".");
    let idx: Vec<&[u8]> = gt_parts(gt).collect();
    if idx.iter().any(|i| matches!(*i, b"." | b"")) {
        return Ok(Call::No(Why::NoCallInGvcf));
    }
    // alleles = [ref] + alts, indexed as a Python list (negative from the end).
    let n_alleles = 1 + r.alts.len() as i64;
    let mut chosen = Vec::with_capacity(idx.len());
    for i in &idx {
        let Some(mut k) = int(i) else {
            return Ok(Call::No(Why::Unreadable));
        };
        if k < 0 {
            k += n_alleles;
        }
        if !(0..n_alleles).contains(&k) {
            return Ok(Call::No(Why::Unreadable));
        }
        chosen.push(if k == 0 {
            r.ref_
        } else {
            r.alts[k as usize - 1]
        });
    }
    if r.pos != pos {
        // The site lies inside a record that starts before it: a reference
        // block (each base is the reference) or a deletion.
        if idx.iter().all(|i| *i == b"0") && r.is_block() {
            return bases(&vec![ref38; idx.len()]);
        }
        return Ok(Call::No(Why::InsideDeletion));
    }
    let Some(&first) = r.ref_.first() else {
        return Err("a record with an empty REF".into());
    };
    if first.to_ascii_uppercase() != ref38 {
        tally.reference_mismatch += 1;
        return Ok(Call::No(Why::RefMismatch));
    }
    tally.reference_checked += 1;
    if !matches!(r.filter, b"PASS" | b".") && idx.iter().any(|i| *i != b"0") {
        return Ok(Call::No(Why::Filtered));
    }
    let mut out = Vec::with_capacity(chosen.len());
    for a in chosen {
        match a.first().map(u8::to_ascii_uppercase) {
            Some(b @ (b'A' | b'C' | b'G' | b'T')) => out.push(b),
            _ => return Ok(Call::No(Why::Symbolic)),
        }
    }
    bases(&out)
}

fn bases(b: &[u8]) -> Result<Call, String> {
    Ok(match *b {
        [a] => Call::Bases { n: 1, b: [a, 0] },
        [a, c] => Call::Bases { n: 2, b: [a, c] },
        _ => Call::Many(u8::try_from(b.len()).unwrap_or(u8::MAX)),
    })
}

/// The walk over one contig's sites.
#[derive(Default)]
struct Walk {
    /// Indices into the table's sites, sorted by GRCh38 position.
    sites: Vec<usize>,
    /// The next site to resolve.
    i: usize,
    /// The lines of the records that start at the last position seen.
    pending: Vec<u8>,
    ends: Vec<usize>,
    /// (first, last) base inside a non-reference longer allele.
    spans: Vec<(i64, i64)>,
    /// The length of `spans` at which the spans behind the walk are dropped.
    prune_at: usize,
}

impl Walk {
    fn push_pending(&mut self, line: &[u8]) {
        self.pending.extend_from_slice(line);
        self.ends.push(self.pending.len());
    }

    fn clear_pending(&mut self) {
        self.pending.clear();
        self.ends.clear();
    }

    /// `resolve(upto)`: gives each site before `upto` (None: every one) its
    /// call from the pending records.
    fn resolve(
        &mut self,
        upto: Option<i64>,
        sites: &mut [Site],
        tally: &mut Tally,
    ) -> Result<(), String> {
        let mut recs = Vec::with_capacity(self.ends.len());
        let mut from = 0;
        for &end in &self.ends {
            if let Some(r) = parse_record(&self.pending[from..end])? {
                recs.push(r);
            }
            from = end;
        }
        let start = recs.first().map(|r| r.pos);
        while self.i < self.sites.len() {
            let s = &mut sites[self.sites[self.i]];
            let p = i64::from(s.pos38);
            if upto.is_some_and(|u| p >= u) {
                break;
            }
            self.i += 1;
            self.spans.retain(|sp| sp.1 >= p);
            if self.spans.iter().any(|&(a, b)| a <= p && p <= b) {
                s.call = Call::No(Why::InsideDeletion);
                continue;
            }
            let hit: Vec<&Rec> = recs
                .iter()
                .filter(|r| start.is_some_and(|st| st <= p) && p <= r.end)
                .collect();
            let Some(&first) = hit.first() else {
                s.call = Call::No(Why::NoRecord);
                continue;
            };
            // At the record's start a variant wins over a block; past it, a block.
            let at_start = Some(p) == start;
            let r = hit
                .iter()
                .find(|r| r.is_block() == at_start)
                .copied()
                .unwrap_or(first);
            s.call = call_from_record(r, p, s.ref38, tally)?;
            if matches!(s.call, Call::Bases { .. } | Call::Many(_)) {
                if r.is_block() || r.pos != p {
                    tally.from_blocks += 1;
                } else {
                    tally.from_variants += 1;
                }
            }
        }
        Ok(())
    }
}

/// The contig of a gVCF line as the site lists are keyed.
fn contig_chrom(contig: &[u8]) -> Option<u8> {
    let c = contig.strip_prefix(b"chr")?;
    if c == b"M" {
        return Some(MT);
    }
    if c == b"MT" {
        return None;
    }
    chrom_index(std::str::from_utf8(c).ok()?)
}

/// Walks the gVCF once and gives each site its call.
///
/// A site's call comes from the records that start last at or before it: a
/// record starting at the site (a variant before a block), else a reference
/// block that covers it. A site inside the span of a variant whose genotype
/// carries a longer allele (a deletion, an MNV) is `--` whatever blocks
/// cover it: DeepVariant starts its next block one base after the variant,
/// inside the deleted bases. Only lines near a site are parsed.
fn genotype_all(
    gvcf: &str,
    input: &mut dyn BufRead,
    table: &mut Table,
    tally: &mut Tally,
) -> Result<(), Fail> {
    let mut by_chrom: Vec<Vec<usize>> = vec![Vec::new(); CHROMS.len()];
    for (k, s) in table.sites.iter().enumerate() {
        by_chrom[usize::from(s.chrom)].push(k);
    }
    for l in &mut by_chrom {
        l.sort_by_key(|&k| table.sites[k].pos38);
    }
    let sites = &mut table.sites;
    let bad = |e: String| Fail::data(format!("{gvcf}: {e}"));
    let mut done: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut cur: Option<Vec<u8>> = None;
    let mut chrom: Option<u8> = None;
    let mut w = Walk::default();
    let mut last_pos: i64 = 0;
    let mut line = Vec::with_capacity(1 << 12);
    let mut lineno = 0u64;
    loop {
        line.clear();
        if input
            .read_until(b'\n', &mut line)
            .map_err(|e| Fail::io(gvcf, &e))?
            == 0
        {
            break;
        }
        lineno += 1;
        // Python's universal newlines: \r\n reads as \n.
        if line.ends_with(b"\r\n") {
            line.truncate(line.len() - 2);
        } else if line.ends_with(b"\n") {
            line.pop();
        }
        if line.first() == Some(&b'#') {
            continue;
        }
        let l = line.as_slice();
        // Python indexed the line with its "\n"; t1 and t2 never reach it.
        let t1 = find(l, b'\t', 0);
        let t2 = find(l, b'\t', t1 + 1);
        let full = l.len() as isize + 1;
        let contig = slice(l, 0, if t1 < 0 { full - 1 } else { t1 });
        let pos_text = slice(l, t1 + 1, if t2 < 0 { full - 1 } else { t2 });
        let pos = int(pos_text)
            .ok_or_else(|| Fail::data(format!("{gvcf} line {lineno}: not a gVCF record")))?;
        if cur.as_deref() != Some(contig) {
            if let Some(old) = cur.take() {
                if let Some(c) = chrom {
                    w.resolve(None, sites, tally).map_err(bad)?;
                    by_chrom[usize::from(c)] = std::mem::take(&mut w.sites);
                }
                done.insert(old);
            }
            if done.contains(contig) {
                return Err(Fail::data(format!(
                    "{gvcf} is not sorted: {} appears twice",
                    String::from_utf8_lossy(contig)
                )));
            }
            cur = Some(contig.to_vec());
            chrom = contig_chrom(contig);
            w = Walk {
                sites: chrom
                    .map(|c| std::mem::take(&mut by_chrom[usize::from(c)]))
                    .unwrap_or_default(),
                ..Walk::default()
            };
            last_pos = 0;
        }
        if w.i >= w.sites.len() {
            continue;
        }
        if pos < last_pos {
            return Err(Fail::data(format!(
                "{gvcf} is not sorted: {}:{pos} after {last_pos}",
                String::from_utf8_lossy(contig)
            )));
        }
        if pos != last_pos {
            // Records start in order, so a site before this one gets its
            // call from the pending records; only then are they parsed.
            if i64::from(sites[w.sites[w.i]].pos38) < pos {
                w.resolve(Some(pos), sites, tally).map_err(bad)?;
            }
            w.clear_pending();
            last_pos = pos;
        }
        w.push_pending(l);
        let t3 = find(l, b'\t', t2 + 1);
        let t4 = find(l, b'\t', t3 + 1);
        if t4 - t3 > 2 {
            let r = parse_record(l).map_err(bad)?;
            if let Some(r) = r
                && !r.is_block()
                && nonref_gt(&r)
            {
                // Every site left starts at or after `pos`, so a span that
                // ends before it covers none: drop those now and then, or a
                // chromosome with few sites keeps every deletion it holds.
                if w.spans.len() >= w.prune_at {
                    w.spans.retain(|sp| sp.1 >= pos);
                    w.prune_at = (2 * w.spans.len()).max(1024);
                }
                w.spans.push((pos + 1, r.end));
            }
        }
    }
    if cur.is_some() {
        w.resolve(None, sites, tally).map_err(bad)?;
    }
    for s in sites.iter_mut() {
        if s.call == Call::Unset {
            s.call = Call::No(Why::NoRecord);
        }
    }
    Ok(())
}

fn complement(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'T' => Some(b'A'),
        _ => None,
    }
}

/// `finalise`: the ploidy rules of the chip files and the GRCh37 strand;
/// each call becomes the letters written, sorted.
fn finalise(table: &mut Table, sex: &str, par: &Par) -> Result<(), Fail> {
    for k in 0..table.sites.len() {
        let s = table.sites[k];
        let mut call = s.call;
        if matches!(call, Call::Bases { .. } | Call::Many(_)) && s.chrom == Y && sex == "female" {
            call = Call::No(Why::ChrYFemale);
        }
        if let Call::Many(n) = call {
            let mut id = String::new();
            table.rsid(s.rsid, &mut id);
            return Err(Fail::data(format!(
                "{id} at {}:{}: a genotype of {n} alleles; a chip file holds one or two",
                gvcf_contig(s.chrom),
                s.pos38
            )));
        }
        if let Call::Bases { n, b } = call {
            let haploid = s.chrom == MT
                || (sex == "male"
                    && (s.chrom == X || s.chrom == Y)
                    && !in_par(par, &gvcf_contig(s.chrom), s.pos38));
            let (mut n, mut b) = (n, b);
            if haploid && n == 2 {
                if b[0] == b[1] {
                    n = 1;
                } else {
                    call = Call::No(Why::HaploidHet);
                }
            }
            if !matches!(call, Call::No(_)) {
                if s.minus {
                    for x in &mut b[..usize::from(n)] {
                        *x = complement(*x).ok_or_else(|| {
                            Fail::data(format!(
                                "the site table's GRCh38 base {} has no complement",
                                py_repr(&char::from(*x).to_string())
                            ))
                        })?;
                    }
                }
                if n == 2 && b[1] < b[0] {
                    b.swap(0, 1);
                }
                call = Call::Bases { n, b };
            }
        }
        table.sites[k].call = call;
    }
    Ok(())
}

/// The letters of a final call, or None for `--`.
fn letters(c: Call) -> Option<([u8; 2], usize)> {
    match c {
        Call::Bases { n, b } => Some((b, usize::from(n))),
        _ => None,
    }
}

/// The AncestryDNA code of a site: 25 for X inside a GRCh37 PAR.
fn ancestry_code(s: &Site) -> u8 {
    if s.chrom == X && PAR37.iter().any(|&(a, b)| a <= s.pos37 && s.pos37 <= b) {
        return 25;
    }
    match s.chrom {
        X => 23,
        Y => 24,
        MT => 26,
        c => c + 1,
    }
}

fn write_23andme(t: &Table, today: &str) -> Vec<u8> {
    let template = t.meta("template_sites").unwrap_or("?");
    let mut out = String::with_capacity(t.sites.len() * 32 + 1024);
    out.push_str("# This data file was generated by OwnGenome (step 40) in the layout of a 23andMe raw data file.\n");
    out.push_str("# It was not made by 23andMe: the genotypes come from whole-genome sequencing calls, not from a chip.\n");
    out.push_str(&format!(
        "# Written {today}. Sites: the 23andMe v5 chip list, {} of {template} sites.\n",
        t.sites.len()
    ));
    out.push_str("# We are using reference human assembly build 37 (also known as GRCh37); genotypes are on its + strand.\n");
    out.push_str("# A genotype of -- means the calls hold no genotype for that site.\n");
    out.push_str(&format!("# More: {DOC_URL}\n"));
    out.push_str("#\n");
    out.push_str(check::HEADER_23ANDME);
    out.push('\n');
    let mut order: Vec<usize> = (0..t.sites.len()).collect();
    order.sort_by(|&a, &b| {
        let (x, y) = (&t.sites[a], &t.sites[b]);
        (x.chrom, x.pos37)
            .cmp(&(y.chrom, y.pos37))
            .then_with(|| t.rsid_cmp(x.rsid, y.rsid))
    });
    let mut id = String::new();
    for k in order {
        let s = &t.sites[k];
        t.rsid(s.rsid, &mut id);
        out.push_str(&id);
        out.push('\t');
        out.push_str(CHROMS[usize::from(s.chrom)]);
        out.push('\t');
        out.push_str(&s.pos37.to_string());
        out.push('\t');
        match letters(s.call) {
            Some((b, n)) => b[..n].iter().for_each(|&c| out.push(char::from(c))),
            None => out.push_str("--"),
        }
        out.push('\n');
    }
    out.into_bytes()
}

fn write_ancestry(t: &Table, today: &str) -> Vec<u8> {
    let template = t.meta("template_sites").unwrap_or("?");
    let mut out = String::with_capacity(t.sites.len() * 34 + 1024);
    out.push_str(check::ANCESTRY_FIRST);
    out.push('\n');
    out.push_str("#This file was generated by OwnGenome (step 40) in the layout of an AncestryDNA raw data file.\n");
    out.push_str("#It was not made by AncestryDNA: the genotypes come from whole-genome sequencing calls, not from a chip,\n");
    out.push_str(&format!(
        "#at the sites of the 23andMe v5 chip list ({} of {template} sites). Written {today}.\n",
        t.sites.len()
    ));
    out.push_str("#Positions are on human assembly build 37 (GRCh37), alleles on its + strand. 0 0 means no genotype.\n");
    out.push_str("#Chromosome codes: 1-22, 23 = X, 24 = Y, 25 = X pseudoautosomal, 26 = MT.\n");
    out.push_str(&format!("#More: {DOC_URL}\n"));
    out.push_str(check::HEADER_ANCESTRY);
    out.push('\n');
    let mut order: Vec<(u8, usize)> = t
        .sites
        .iter()
        .enumerate()
        .map(|(k, s)| (ancestry_code(s), k))
        .collect();
    order.sort_by(|&(ca, a), &(cb, b)| {
        let (x, y) = (&t.sites[a], &t.sites[b]);
        (ca, x.pos37)
            .cmp(&(cb, y.pos37))
            .then_with(|| t.rsid_cmp(x.rsid, y.rsid))
    });
    let mut id = String::new();
    for (code, k) in order {
        let s = &t.sites[k];
        t.rsid(s.rsid, &mut id);
        let (a, b) = match letters(s.call) {
            Some((b, n)) => (char::from(b[0]), char::from(b[n - 1])),
            None => ('0', '0'),
        };
        out.push_str(&format!("{id}\t{code}\t{}\t{a}\t{b}\n", s.pos37));
    }
    out.into_bytes()
}

/// The date the header carries and the zip members are stamped with.
pub fn date_of(flag: Option<&str>, epoch: Option<&str>, now: u64) -> Result<DosDate, String> {
    if let Some(d) = flag {
        return DosDate::parse(d).ok_or_else(|| {
            format!(
                "argument --date: {} is not a date from 1980-01-01 to 2107-12-31 (YYYY-MM-DD)",
                py_repr(d)
            )
        });
    }
    let secs = match epoch {
        Some(e) => e.trim().parse::<u64>().map_err(|_| {
            format!(
                "SOURCE_DATE_EPOCH={} is not a number of seconds",
                py_repr(e)
            )
        })?,
        None => now,
    };
    DosDate::from_unix_days(secs / 86_400)
        .ok_or_else(|| "the date is outside 1980-01-01 to 2107-12-31".to_owned())
}

/// The options of `export`.
pub struct ExportArgs<'a> {
    pub gvcf: &'a str,
    pub sites: &'a str,
    pub sample: &'a str,
    pub sex: &'a str,
    pub par: Option<&'a str>,
    pub out_dir: &'a str,
    pub date: DosDate,
}

/// Runs the export; the stdout summary on success.
pub fn export(a: &ExportArgs, workers: NonZero<usize>) -> Result<String, Fail> {
    let mut table = read_table(a.sites)?;
    let par = read_par(a.par)?;
    let mut tally = Tally::default();
    let mut input = og_io::bgzf::open(a.gvcf, workers).map_err(|e| Fail::io(a.gvcf, &e))?;
    genotype_all(a.gvcf, &mut input, &mut table, &mut tally)?;
    drop(input);
    finalise(&mut table, a.sex, &par)?;
    let (checked, bad) = (tally.reference_checked, tally.reference_mismatch);
    if bad > 10.max((checked + bad) / 100) {
        return Err(Fail::data(format!(
            "{bad} of {} gVCF records at a site have another reference base than the site table: \
             the gVCF is not on GRCh38, or the site table is wrong. Nothing was written.",
            checked + bad
        )));
    }
    std::fs::create_dir_all(a.out_dir).map_err(|e| Fail::io(a.out_dir, &e))?;
    let today = a.date.to_string();
    let base = path_join(a.out_dir, a.sample);
    let (t23, tanc) = (
        format!("{base}_23andme_v5.txt"),
        format!("{base}_ancestrydna_layout.txt"),
    );
    // One text at a time: built, read back, zipped in memory, dropped.
    type Layouts<'a> = [(&'a str, fn(&Table, &str) -> Vec<u8>, Layout, String); 2];
    let layouts: Layouts = [
        (
            &t23,
            write_23andme,
            Layout::TwentyThree,
            format!("{base}_23andme_v5.zip"),
        ),
        (
            &tanc,
            write_ancestry,
            Layout::Ancestry,
            format!("{base}_ancestrydna_layout.zip"),
        ),
    ];
    let mut problems = Vec::new();
    let mut zips = Vec::new();
    for (txt, write, layout, zpath) in &layouts {
        let text = write(&table, &today);
        problems.extend(check::check_text(txt, &text, *layout));
        if problems.is_empty() {
            let mut z = Vec::new();
            zip::write_one(&mut z, basename(txt), &text, a.date)
                .map_err(|e| Fail::io(zpath, &e))?;
            zips.push((zpath, z));
        }
    }
    if !problems.is_empty() {
        for (txt, write, _, _) in &layouts {
            std::fs::write(txt, write(&table, &today)).map_err(|e| Fail::io(txt, &e))?;
        }
        for p in problems.iter().take(20) {
            eprintln!("  {p}");
        }
        return Err(Fail::data(format!(
            "the written files fail the strict reader ({} problems); they are left unzipped in {}",
            problems.len(),
            a.out_dir
        )));
    }
    for (zpath, z) in zips {
        let mut f = AtomicFile::create(zpath).map_err(|e| Fail::io(zpath, &e))?;
        f.write_all(&z)
            .and_then(|()| f.commit())
            .map_err(|e| Fail::io(zpath, &e))?;
    }

    let called = table
        .sites
        .iter()
        .filter(|s| matches!(s.call, Call::Bases { .. }))
        .count();
    let n = table.sites.len();
    let mut reasons: Vec<(String, u64)> = Vec::new();
    for s in &table.sites {
        if let Call::No(why) = s.call {
            let key = format!("no_call: {why}");
            match reasons.iter_mut().find(|(k, _)| *k == key) {
                Some(r) => r.1 += 1,
                None => reasons.push((key, 1)),
            }
        }
    }
    reasons.sort();
    let stem = basename(&base);
    let mut tsv = String::from("key\tvalue\n");
    let rows: [(&str, String); 12] = [
        ("sites", n.to_string()),
        ("called", called.to_string()),
        ("no_call", (n - called).to_string()),
        ("from_reference_blocks", tally.from_blocks.to_string()),
        ("from_variant_records", tally.from_variants.to_string()),
        ("reference_checked", checked.to_string()),
        ("reference_mismatch", bad.to_string()),
        (
            "template_sites",
            table.meta("template_sites").unwrap_or(".").to_owned(),
        ),
        (
            "table_dropped",
            table.meta("dropped").unwrap_or(".").to_owned(),
        ),
        (
            "sex",
            if a.sex.is_empty() { "not given" } else { a.sex }.to_owned(),
        ),
        ("file_23andme", format!("{stem}_23andme_v5.zip")),
        ("file_ancestrydna", format!("{stem}_ancestrydna_layout.zip")),
    ];
    for (k, v) in &rows {
        tsv.push_str(&format!("{k}\t{v}\n"));
    }
    for (k, v) in &reasons {
        tsv.push_str(&format!("{k}\t{v}\n"));
    }
    let summary = format!("{base}_consumer_export.tsv");
    let mut f = AtomicFile::create(&summary).map_err(|e| Fail::io(&summary, &e))?;
    f.write_all(tsv.as_bytes())
        .and_then(|()| f.commit())
        .map_err(|e| Fail::io(&summary, &e))?;

    let pct = 100.0 * called as f64 / n as f64;
    let mut out = format!(
        "Sites: {n} (of {} on the chip list; {} dropped when the table was built)\n",
        table.meta("template_sites").unwrap_or("?"),
        table.meta("dropped").unwrap_or("?")
    );
    out.push_str(&format!(
        "Called: {called} ({}%), {} of them from reference blocks; no-call (--): {}\n",
        pycompat::fmt::fixed(pct, 1),
        tally.from_blocks,
        n - called
    ));
    for (k, v) in &reasons {
        out.push_str(&format!("  {}: {v}\n", &k["no_call: ".len()..]));
    }
    out.push_str(&format!(
        "Reference bases checked at {checked} variant records, {bad} differ from the table\n"
    ));
    Ok(out)
}

/// `og consumer export ...`.
pub fn cli(args: &[String]) -> Result<u8, Usage> {
    let p = Args::parse(
        args,
        &["gvcf", "sites", "sample", "sex", "par", "out-dir", "date"],
    )?;
    if let Some(extra) = p.positional.first() {
        return Err(Usage(format!("unrecognized arguments: {extra}")));
    }
    let sex = p.get("sex").unwrap_or("");
    if !matches!(sex, "" | "male" | "female") {
        return Err(Usage(format!(
            "argument --sex: invalid choice: {} (choose from 'male', 'female', '')",
            py_repr(sex)
        )));
    }
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    let epoch = std::env::var("SOURCE_DATE_EPOCH").ok();
    let date = date_of(p.get("date"), epoch.as_deref(), now)?;
    let a = ExportArgs {
        gvcf: p.need("gvcf")?,
        sites: p.need("sites")?,
        sample: p.need("sample")?,
        sex,
        par: p.get("par"),
        out_dir: p.need("out-dir")?,
        date,
    };
    let workers = og_io::bgzf::workers().min(NonZero::new(4).expect("4 > 0"));
    Ok(match export(&a, workers) {
        Ok(summary) => {
            print!("{summary}");
            0
        }
        Err(f) => report(&f),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// tests/parity/inputs/consumer: the cases of the old unit test.
    fn input(name: &str) -> String {
        let p: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../../../tests/parity/inputs/consumer",
            name,
        ]
        .iter()
        .collect();
        p.to_string_lossy().into_owned()
    }

    /// Runs the export into a fresh directory; the 23andMe calls keyed by
    /// rsID, and the summary table.
    fn run(gvcf: &str, sites: &str, sex: &str) -> Result<(Vec<(String, String)>, String), Fail> {
        let out = std::env::temp_dir().join(format!(
            "og-consumer-test-{}-{gvcf}-{sites}-{sex}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&out);
        let out_dir = out.to_string_lossy().into_owned();
        let a = ExportArgs {
            gvcf: &input(gvcf),
            sites: &input(sites),
            sample: "S",
            sex,
            par: Some(&input("par_grch38.bed")),
            out_dir: &out_dir,
            date: DosDate::parse("2026-01-02").expect("a date"),
        };
        let res = export(&a, NonZero::new(2).expect("2 > 0"));
        let read = || {
            let z = std::fs::read(out.join("S_23andme_v5.zip")).expect("the zip");
            let members = zip::read_all(&z).expect("a zip");
            let text = String::from_utf8(members[0].1.clone()).expect("ASCII");
            let calls = text
                .lines()
                .filter(|l| !l.starts_with('#'))
                .map(|l| {
                    let c: Vec<&str> = l.split('\t').collect();
                    (c[0].to_owned(), c[3].to_owned())
                })
                .collect();
            let tsv = std::fs::read_to_string(out.join("S_consumer_export.tsv")).expect("tsv");
            (calls, tsv)
        };
        let got = res.map(|_| read());
        let _ = std::fs::remove_dir_all(&out);
        got
    }

    fn call<'a>(calls: &'a [(String, String)], id: &str) -> &'a str {
        &calls.iter().find(|(r, _)| r == id).expect("the site").1
    }

    #[test]
    fn every_site_has_the_genotype_the_gvcf_holds() {
        // (site number, genotype) of a male: a block is reference, a het and
        // a hom-alt call keep their bases, a filtered call, a no-call block,
        // a gap, a site inside a called deletion (also under the block that
        // starts inside it) and a contig the gVCF lacks are '--'; a block
        // over a reference-called deletion is reference; a variant wins
        // over a block at its position; strand '-' is complemented; chrX
        // outside the PAR gets one letter ('--' when het), inside two;
        // chrY and chrM one.
        let want = [
            (101, "AA"),
            (102, "CT"),
            (103, "AA"),
            (104, "--"),
            (105, "--"),
            (106, "--"),
            (107, "--"),
            (108, "AA"),
            (109, "CG"),
            (110, "--"),
            (111, "CC"),
            (201, "AG"),
            (202, "T"),
            (203, "--"),
            (204, "G"),
            (301, "C"),
            (401, "A"),
            (402, "--"),
        ];
        let (calls, tsv) = run("unit.g.vcf.gz", "unit-sites.tsv.gz", "male").expect("export");
        for (n, g) in want {
            assert_eq!(call(&calls, &format!("rs{n}")), g, "rs{n}");
        }
        for row in [
            "sites\t18",
            "called\t11",
            "no_call\t7",
            "no_call: filtered call\t1",
            "no_call: no gVCF record\t2",
            "no_call: no-call in the gVCF\t1",
            "no_call: inside a deletion or a longer allele\t1",
            "reference_mismatch\t0",
        ] {
            assert!(tsv.lines().any(|l| l == row), "{row} in {tsv}");
        }
        let (female, _) = run("unit.g.vcf.gz", "unit-sites.tsv.gz", "female").expect("export");
        // Ids and genotypes apart, so no line pairs an rsID with letters.
        let ids = ["rs301", "rs202", "rs203", "rs204"];
        let want = ["--", "TT", "CT", "G"];
        for (id, g) in ids.iter().zip(want) {
            assert_eq!(call(&female, id), g, "female {id}");
        }
        let (holed, _) = run("hole.g.vcf.gz", "unit-sites.tsv.gz", "male").expect("export");
        assert_eq!(call(&holed, "rs101"), "--");
        assert_eq!(
            holed
                .iter()
                .filter(|(id, g)| id != "rs101" && call(&calls, id) != g)
                .count(),
            0
        );
    }

    #[test]
    fn a_table_off_by_one_or_an_unsorted_gvcf_stops_the_export() {
        let (_, tsv) = run("run.g.vcf.gz", "run-sites.tsv.gz", "").expect("export");
        assert!(tsv.lines().any(|l| l == "called\t50"), "{tsv}");
        let e = run("run.g.vcf.gz", "shifted-sites.tsv.gz", "").unwrap_err();
        assert!(
            e.code == 1 && e.msg.contains("another reference base than the site table"),
            "{e:?}"
        );
        let e = run("unsorted.g.vcf.gz", "unit-sites.tsv.gz", "").unwrap_err();
        assert!(e.code == 1 && e.msg.contains("is not sorted"), "{e:?}");
    }

    #[test]
    fn rsids_keep_their_spelling() {
        let mut names = Vec::new();
        assert_eq!(RsId::new("rs10", &mut names), RsId::Rs(10));
        assert_eq!(RsId::new("i7001", &mut names), RsId::I(7001));
        assert_eq!(RsId::new("rs0", &mut names), RsId::Rs(0));
        assert_eq!(RsId::new("rs010", &mut names), RsId::Name(0));
        assert_eq!(RsId::new("rs99999999999", &mut names), RsId::Name(1));
        assert_eq!(&*names[0], "rs010");
    }

    #[test]
    fn python_find_and_slice() {
        assert_eq!(find(b"a\tb", b'\t', 0), 1);
        assert_eq!(find(b"a\tb", b'\t', 2), -1);
        assert_eq!(find(b"ab", b'\t', 5), -1);
        assert_eq!(slice(b"abc", 0, -1), b"ab");
        assert_eq!(slice(b"abc", 2, 1), b"");
    }

    #[test]
    fn end_as_re_search_finds_it() {
        assert_eq!(info_end(b"END=120"), Some(120));
        assert_eq!(info_end(b"X=1;END=5;Y"), Some(5));
        assert_eq!(info_end(b"XEND=9;END=7"), Some(7));
        assert_eq!(info_end(b"END=;END=3"), Some(3));
        assert_eq!(info_end(b"."), None);
    }

    #[test]
    fn dates() {
        assert_eq!(
            date_of(Some("2026-11-10"), None, 0).unwrap().to_string(),
            "2026-11-10"
        );
        assert!(date_of(Some("2026-13-01"), None, 0).is_err());
        assert_eq!(
            date_of(None, Some("0"), 9).map(|d| d.to_string()),
            Err("the date is outside 1980-01-01 to 2107-12-31".to_owned())
        );
        assert_eq!(
            date_of(None, Some("1762732800"), 0).unwrap().to_string(),
            "2025-11-10"
        );
        assert_eq!(
            date_of(None, None, 1_762_819_199).unwrap().to_string(),
            "2025-11-10"
        );
    }
}
