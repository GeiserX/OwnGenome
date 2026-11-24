//! `og sv consensus-check`: each record of step 22's SURVIVOR merge checked
//! so that a consensus record is one event, called by two or more callers.
//!
//! SURVIVOR 1.0.7 (`merge LIST 1000 2 1 1 0 50`) pairs calls whose
//! breakpoints lie within 1,000 bp of each other. That alone merges two
//! deletions of 100 bp that start 900 bp apart and do not overlap, and on
//! some inputs it pairs a deletion with a duplication at the same place
//! although the type and strand flags are set. The record then reads "2
//! callers" for an event no two callers agree on.
//!
//! The members of a record are read from its caller columns: FORMAT `CO`
//! (each member's coordinates, comma-separated when one caller gave two
//! calls), `TY` (types), `ST` (strands) and `LN` (length). Two members are
//! the same event when all of these hold:
//!
//! - the same type, by SURVIVOR's own rule: `INV`, `BND` and `TRA` count as
//!   one type, as they do in SURVIVOR; any other type only matches itself
//!   (`DUP:TANDEM` is a `DUP`);
//! - the same strands (`ST`);
//! - insertions: starts within [`INS_MAX_DIST`] bp, and the shorter at
//!   least [`INS_MIN_RATIO`] of the longer (`LN`);
//! - a breakend to another contig: the same contig, both breakpoints
//!   within [`BND_MAX_DIST`] bp (SURVIVOR's own rule, as there is no span);
//! - anything else: a reciprocal overlap of at least [`MIN_OVERLAP`], that
//!   is, an overlap of at least that share of the longer call.
//!
//! The members are grouped in column order (the order of SUPP_VEC): each
//! joins the first group whose every member it matches, or starts a new
//! one. A record whose members form one group is written as SURVIVOR wrote
//! it. Otherwise each group with two or more callers becomes a record of
//! its own and the rest are dropped. Such a record gets SUPP, SUPP_VEC,
//! SVTYPE, STRANDS and CHR2 of its members, POS and END of its middle
//! member (by start; the lower of the two middle ones), SVLEN as the mean
//! of the members' lengths (negative for a deletion, SURVIVOR's own way),
//! CIPOS and CIEND as the spread of the members' starts and ends around
//! POS and END, REF `N` and a symbolic ALT. A caller column without a
//! member gets SURVIVOR's empty column; a column that keeps only some of
//! its members keeps those entries of `CO`, `TY`, `ST` and `QV`, and its
//! `ID` only when it names a kept member (SURVIVOR writes the ID of the
//! last one), `.` otherwise.
//!
//! A record whose columns cannot be read this way (no `CO`, a coordinate
//! that does not parse, a SUPP_VEC that does not fit the columns) is
//! written as SURVIVOR wrote it and counted apart. A summary line goes to
//! stderr.

use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};

use og_filters::args::{self, ArgError, Spec};

use super::info_value;

const USAGE: &str = "usage: og sv consensus-check IN.vcf[.gz]|- OUT.vcf|-
";

const SPEC: Spec = Spec {
    prog: "og sv consensus-check",
    usage: USAGE,
    options: &[],
    positionals: &["IN", "OUT"],
    required: &[],
};

/// The reciprocal overlap two calls of a span need (DEL, DUP, INV).
pub const MIN_OVERLAP: f64 = 0.5;
/// How far apart two insertions may start, in bp.
pub const INS_MAX_DIST: i64 = 1000;
/// The share of the longer insertion the shorter must reach.
pub const INS_MIN_RATIO: f64 = 0.5;
/// How far apart each breakpoint of two breakends to another contig may
/// be, in bp: SURVIVOR's own distance in step 22.
pub const BND_MAX_DIST: i64 = 1000;

/// The header line the check adds, just above `#CHROM`.
pub const HEADER: &str = "##og_sv_consensus_check=<Description=\"Each record is one event two or more callers agree on: same type (INV, BND and TRA as one), same strands, reciprocal overlap >= 0.5; insertions start within 1000 bp, the shorter >= 0.5 of the longer; breakends to another contig within 1000 bp at both ends. Records that failed were split by event, and dropped when fewer than 2 callers were left.\">";

/// SURVIVOR's column for a caller with no member.
const EMPTY: &str = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
const SURVIVOR_FORMAT: &str = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO";

/// One call of one caller inside a merged record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Member {
    /// The caller's column, 0 for the first.
    pub caller: usize,
    /// Its place in the column's `CO` list.
    pub index: usize,
    /// `TY` (`DEL`, `DUP`, `INV`, `TRA`, `INS`, ...).
    pub kind: String,
    /// `ST` (`+-`, `-+`, `++`, `--`).
    pub strands: String,
    /// The contig of the second breakpoint.
    pub chrom2: String,
    /// Start and end on the record's contig; for a breakend to another
    /// contig, the position on each.
    pub start: i64,
    pub end: i64,
    /// The length: `LN` for an insertion, 0 for a breakend to another
    /// contig, the span otherwise.
    pub len: i64,
    /// Whether the second breakpoint is on another contig.
    pub other_contig: bool,
}

/// Why two members are not one event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Apart {
    /// Different types, by SURVIVOR's rule.
    Type,
    /// Different strands.
    Strands,
    /// Spans that overlap by less than [`MIN_OVERLAP`].
    Overlap,
    /// Insertions more than [`INS_MAX_DIST`] apart, or breakends to
    /// another contig more than [`BND_MAX_DIST`] apart or on two contigs.
    Distance,
    /// Insertions whose sizes differ more than [`INS_MIN_RATIO`] allows.
    Size,
}

/// A type as SURVIVOR's rule compares it: the part before `:`, with
/// `INV`, `BND` and `TRA` as one.
fn type_class(kind: &str) -> &str {
    let base = kind.split(':').next().unwrap_or(kind);
    match base {
        "INV" | "BND" | "TRA" => "INV",
        other => other,
    }
}

/// Whether two members are one event (`None`), or why not.
pub fn apart(a: &Member, b: &Member) -> Option<Apart> {
    if type_class(&a.kind) != type_class(&b.kind) {
        return Some(Apart::Type);
    }
    if a.strands != b.strands {
        return Some(Apart::Strands);
    }
    if type_class(&a.kind) == "INS" {
        if (a.start - b.start).abs() > INS_MAX_DIST {
            return Some(Apart::Distance);
        }
        let (x, y) = (a.len.abs(), b.len.abs());
        let (short, long) = (x.min(y), x.max(y));
        if long == 0 || (short as f64) < INS_MIN_RATIO * long as f64 {
            return Some(Apart::Size);
        }
        return None;
    }
    if a.other_contig || b.other_contig {
        if a.chrom2 != b.chrom2
            || (a.start - b.start).abs() > BND_MAX_DIST
            || (a.end - b.end).abs() > BND_MAX_DIST
        {
            return Some(Apart::Distance);
        }
        return None;
    }
    let overlap = a.end.min(b.end) - a.start.max(b.start);
    let longer = (a.end - a.start).max(b.end - b.start).max(1);
    if overlap <= 0 || (overlap as f64) < MIN_OVERLAP * longer as f64 {
        return Some(Apart::Overlap);
    }
    None
}

/// The members grouped into events: each member joins the first group
/// whose every member it matches, or starts a new one. Indexes into
/// `members`, in their order.
pub fn group(members: &[Member]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (i, m) in members.iter().enumerate() {
        match groups
            .iter_mut()
            .find(|g| g.iter().all(|&j| apart(&members[j], m).is_none()))
        {
            Some(g) => g.push(i),
            None => groups.push(vec![i]),
        }
    }
    groups
}

/// `CHROM_POS` of a SURVIVOR coordinate, split at the last `_`.
fn contig_pos(s: &str) -> Option<(&str, i64)> {
    let (chrom, pos) = s.rsplit_once('_')?;
    if chrom.is_empty() || pos.is_empty() || !pos.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    Some((chrom, pos.parse().ok()?))
}

/// One `CO` entry, `CHROM_POS-CHROM2_POS2`, on the record's contig
/// `chrom` (contig names may hold `_` and `-`).
fn coords(entry: &str, chrom: &str) -> Option<(String, i64, i64)> {
    let rest = entry.strip_prefix(chrom)?.strip_prefix('_')?;
    let dash = rest.find('-')?;
    let (pos1, second) = (&rest[..dash], &rest[dash + 1..]);
    if pos1.is_empty() || !pos1.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let (chrom2, pos2) = contig_pos(second)?;
    Some((chrom2.to_owned(), pos1.parse().ok()?, pos2))
}

/// Whether a column value says the caller has no member.
fn absent(v: &str) -> bool {
    matches!(v, "" | "." | "NAN" | "NaN" | "NA")
}

/// A record split into its columns, with what the check needs.
struct Record<'a> {
    fields: Vec<&'a str>,
    keys: Vec<&'a str>,
    members: Vec<Member>,
    /// Members per column.
    per_column: Vec<usize>,
}

impl<'a> Record<'a> {
    fn key(&self, k: &str) -> Option<usize> {
        self.keys.iter().position(|&x| x == k)
    }

    /// The record's members, or `None` when its columns cannot be read.
    fn read(line: &'a str) -> Option<Record<'a>> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            return None;
        }
        let keys: Vec<&str> = fields[8].split(':').collect();
        let at = |k: &str| keys.iter().position(|&x| x == k);
        let (co, ty, st, ln_at) = (at("CO")?, at("TY")?, at("ST")?, at("LN")?);
        let chrom = fields[0];
        let bits = info_value(fields[7], "SUPP_VEC").as_bytes();
        let columns = &fields[9..];
        if bits.len() != columns.len() {
            return None;
        }
        let mut members = Vec::new();
        let mut per_column = Vec::with_capacity(columns.len());
        for (c, (col, &bit)) in columns.iter().zip(bits).enumerate() {
            let vals: Vec<&str> = col.split(':').collect();
            let get = |i: usize| vals.get(i).copied().unwrap_or(".");
            let entries = get(co);
            if absent(entries) {
                if bit == b'1' {
                    return None;
                }
                per_column.push(0);
                continue;
            }
            if bit != b'1' {
                return None;
            }
            let kinds: Vec<&str> = get(ty).split(',').collect();
            let strands: Vec<&str> = get(st).split(',').collect();
            let ln: i64 = get(ln_at).parse().map(|v: i64| v.abs()).unwrap_or(-1);
            let mut n = 0;
            for (i, e) in entries.split(',').enumerate() {
                let (chrom2, p1, p2) = coords(e, chrom)?;
                let kind = kinds.get(i).or(kinds.last())?.to_string();
                let strand = strands.get(i).or(strands.last())?.to_string();
                if absent(&kind) || absent(&strand) {
                    return None;
                }
                let other_contig = chrom2 != chrom;
                let (start, end) = if other_contig {
                    (p1, p2)
                } else {
                    (p1.min(p2), p1.max(p2))
                };
                let len = if type_class(&kind) == "INS" {
                    if ln < 0 {
                        return None;
                    }
                    ln
                } else if other_contig {
                    0
                } else {
                    end - start
                };
                members.push(Member {
                    caller: c,
                    index: i,
                    kind,
                    strands: strand,
                    chrom2,
                    start,
                    end,
                    len,
                    other_contig,
                });
                n += 1;
            }
            per_column.push(n);
        }
        if members.is_empty() {
            return None;
        }
        Some(Record {
            fields,
            keys,
            members,
            per_column,
        })
    }

    /// The record of one group of members, as text.
    fn rebuild(&self, g: &[usize]) -> String {
        let kept: Vec<&Member> = g.iter().map(|&i| &self.members[i]).collect();
        let n_cols = self.per_column.len();
        let mut supp_vec = vec![b'0'; n_cols];
        for m in &kept {
            supp_vec[m.caller] = b'1';
        }
        let supp = supp_vec.iter().filter(|&&b| b == b'1').count();
        // The middle member by start; the lower of the two middle ones.
        let mut by_start: Vec<&Member> = kept.clone();
        by_start.sort_by_key(|m| (m.start, m.end));
        let mid = by_start[(by_start.len() - 1) / 2];
        let (pos, end) = (mid.start, mid.end);
        let min_s = kept.iter().map(|m| m.start).min().unwrap_or(pos);
        let max_s = kept.iter().map(|m| m.start).max().unwrap_or(pos);
        let min_e = kept.iter().map(|m| m.end).min().unwrap_or(end);
        let max_e = kept.iter().map(|m| m.end).max().unwrap_or(end);
        let first = kept[0];
        let old_type = info_value(self.fields[7], "SVTYPE");
        // The record's SVTYPE stays when it is the members' type (a BND
        // pair SURVIVOR called INV stays INV); else the first member's.
        let svtype = if !old_type.is_empty() && type_class(old_type) == type_class(&first.kind) {
            old_type.to_owned()
        } else {
            first
                .kind
                .split(':')
                .next()
                .unwrap_or(&first.kind)
                .to_owned()
        };
        let mean = kept.iter().map(|m| m.len.abs()).sum::<i64>() as f64 / kept.len() as f64;
        let mut svlen = mean.round() as i64;
        if type_class(&svtype) == "DEL" {
            svlen = -svlen;
        }
        let info = set_info(
            self.fields[7],
            &[
                ("SUPP", supp.to_string()),
                ("SUPP_VEC", String::from_utf8_lossy(&supp_vec).into_owned()),
                ("SVLEN", svlen.to_string()),
                ("SVTYPE", svtype.clone()),
                ("CHR2", first.chrom2.clone()),
                ("END", end.to_string()),
                ("CIPOS", format!("{},{}", min_s - pos, max_s - pos)),
                ("CIEND", format!("{},{}", min_e - end, max_e - end)),
                ("STRANDS", first.strands.clone()),
            ],
        );
        let mut cols: Vec<String> = Vec::with_capacity(n_cols);
        for (c, col) in self.fields[9..].iter().enumerate() {
            let mine: Vec<usize> = kept
                .iter()
                .filter(|m| m.caller == c)
                .map(|m| m.index)
                .collect();
            if mine.is_empty() {
                cols.push(self.empty_column());
            } else if mine.len() == self.per_column[c] {
                cols.push((*col).to_owned());
            } else {
                cols.push(self.part_column(col, &mine, c));
            }
        }
        let id = self.middle_id(mid);
        let mut out = vec![
            self.fields[0].to_owned(),
            pos.to_string(),
            id,
            "N".to_owned(),
            format!("<{svtype}>"),
            self.fields[5].to_owned(),
            self.fields[6].to_owned(),
            info,
            self.fields[8].to_owned(),
        ];
        out.extend(cols);
        out.join("\t")
    }

    /// The source ID of the middle member, when its column's ID names it
    /// (SURVIVOR writes the last member's), else `.`.
    fn middle_id(&self, mid: &Member) -> String {
        let Some(k) = self.key("ID") else {
            return ".".to_owned();
        };
        if mid.index + 1 != self.per_column[mid.caller] {
            return ".".to_owned();
        }
        match self.fields[9 + mid.caller].split(':').nth(k) {
            Some(v) if !absent(v) => v.to_owned(),
            _ => ".".to_owned(),
        }
    }

    /// SURVIVOR's column for a caller with no member, or `./.` and `.`
    /// for each key of another FORMAT.
    fn empty_column(&self) -> String {
        if self.fields[8] == SURVIVOR_FORMAT {
            return EMPTY.to_owned();
        }
        let mut v = vec!["."; self.keys.len().max(1)];
        if self.keys.first() == Some(&"GT") {
            v[0] = "./.";
        }
        v.join(":")
    }

    /// A column that keeps the members `mine` (indexes into its `CO`) of
    /// its `n` members: those entries of the per-member lists, its ID only
    /// when the last member is kept, LN as the kept members' mean length.
    fn part_column(&self, col: &str, mine: &[usize], c: usize) -> String {
        let n = self.per_column[c];
        let mut vals: Vec<String> = col.split(':').map(str::to_owned).collect();
        for k in ["CO", "TY", "ST", "QV"] {
            if let Some(i) = self.key(k)
                && let Some(v) = vals.get_mut(i)
            {
                let list: Vec<&str> = v.split(',').collect();
                if list.len() == n {
                    *v = mine.iter().map(|&j| list[j]).collect::<Vec<_>>().join(",");
                }
            }
        }
        if let Some(i) = self.key("ID")
            && let Some(v) = vals.get_mut(i)
            && !mine.contains(&(n - 1))
        {
            *v = ".".to_owned();
        }
        if let Some(i) = self.key("LN")
            && let Some(v) = vals.get_mut(i)
        {
            let lens: Vec<i64> = self
                .members
                .iter()
                .filter(|m| m.caller == c && mine.contains(&m.index))
                .map(|m| m.len.abs())
                .collect();
            let mean = lens.iter().sum::<i64>() as f64 / lens.len().max(1) as f64;
            *v = (mean.round() as i64).to_string();
        }
        vals.join(":")
    }
}

/// INFO with the values of `set` in place of the keys it names; keys the
/// INFO does not have are not added.
fn set_info(info: &str, set: &[(&str, String)]) -> String {
    info.split(';')
        .map(|kv| {
            let k = kv.split_once('=').map_or(kv, |(k, _)| k);
            match set.iter().find(|(key, _)| *key == k) {
                Some((key, v)) => format!("{key}={v}"),
                None => kv.to_owned(),
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// What one run read and wrote.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Counts {
    /// Records read.
    pub records: usize,
    /// Records whose members are one event, written as SURVIVOR wrote them.
    pub kept: usize,
    /// Records whose members are not one event.
    pub split: usize,
    /// Records the split ones gave (events of two or more callers).
    pub from_split: usize,
    /// Split records that left no event of two or more callers.
    pub dropped: usize,
    /// Records whose columns could not be read, written as they were.
    pub unread: usize,
}

impl Counts {
    /// Records written.
    pub fn written(&self) -> usize {
        self.kept + self.from_split + self.unread
    }
}

impl fmt::Display for Counts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} records, {} written; {} one event as merged, {} not one event: split into {} records of 2+ callers, {} dropped with fewer than 2 callers left",
            self.records,
            self.written(),
            self.kept,
            self.split,
            self.from_split,
            self.dropped
        )?;
        if self.unread > 0 {
            write!(
                f,
                "; {} could not be read and are written as SURVIVOR wrote them",
                self.unread
            )?;
        }
        Ok(())
    }
}

/// A failure that ends the run with exit 1.
#[derive(Debug)]
pub enum Error {
    Read(String, io::Error),
    Write(String, io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Read(p, e) => write!(f, "cannot read {p}: {e}"),
            Error::Write(p, e) => write!(f, "cannot write {p}: {e}"),
        }
    }
}

/// The records one input record becomes, and how it counts.
pub fn check_line(line: &str, c: &mut Counts) -> Vec<String> {
    c.records += 1;
    let Some(rec) = Record::read(line) else {
        c.unread += 1;
        return vec![line.to_owned()];
    };
    let groups = group(&rec.members);
    let two_callers = |g: &Vec<usize>| {
        let mut callers: Vec<usize> = g.iter().map(|&i| rec.members[i].caller).collect();
        callers.sort_unstable();
        callers.dedup();
        callers.len() >= 2
    };
    if groups.len() == 1 && two_callers(&groups[0]) {
        c.kept += 1;
        return vec![line.to_owned()];
    }
    c.split += 1;
    let out: Vec<String> = groups
        .iter()
        .filter(|g| two_callers(g))
        .map(|g| rec.rebuild(g))
        .collect();
    if out.is_empty() {
        c.dropped += 1;
    }
    c.from_split += out.len();
    out
}

/// Reads SURVIVOR's merge `input` (named `path`, for messages) and writes
/// the checked records to `out`.
pub fn check(input: impl BufRead, path: &str, out: &mut impl Write) -> Result<Counts, Error> {
    let w = |e| Error::Write("the output".to_owned(), e);
    let mut c = Counts::default();
    for line in input.lines() {
        let line = line.map_err(|e| Error::Read(path.to_owned(), e))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        if line.starts_with("#CHROM") {
            writeln!(out, "{HEADER}").map_err(w)?;
        }
        if line.starts_with('#') {
            writeln!(out, "{line}").map_err(w)?;
            continue;
        }
        for rec in check_line(line, &mut c) {
            writeln!(out, "{rec}").map_err(w)?;
        }
    }
    out.flush().map_err(w)?;
    Ok(c)
}

/// Runs the verb on files: `out_path` `-` is stdout.
fn run(input: &str, out_path: &str, stdout: &mut dyn Write) -> Result<Counts, Error> {
    let reader = og_io::open(input).map_err(|e| Error::Read(input.to_owned(), e))?;
    if out_path == "-" {
        let mut w = BufWriter::new(stdout);
        return check(reader, input, &mut w);
    }
    // Written beside the target and renamed, so a failed run leaves no
    // half file behind.
    let part = format!("{out_path}.part");
    let file = File::create(&part).map_err(|e| Error::Write(part.clone(), e))?;
    let mut w = BufWriter::new(file);
    let c = check(reader, input, &mut w).inspect_err(|_| {
        let _ = std::fs::remove_file(&part);
    })?;
    drop(w);
    std::fs::rename(&part, out_path).map_err(|e| Error::Write(out_path.to_owned(), e))?;
    Ok(c)
}

/// `og sv consensus-check IN OUT`: exit 0, 1 for an input or output
/// problem, 2 for a usage error.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let p = match args::parse(&SPEC, argv) {
        Ok(p) => p,
        Err(ArgError::Help) => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        Err(ArgError::Usage(m)) => {
            let _ = writeln!(err, "{USAGE}og sv consensus-check: error: {m}");
            return 2;
        }
    };
    match run(&p.pos[0], &p.pos[1], out) {
        Ok(c) => {
            let _ = writeln!(err, "og sv consensus-check: {c}");
            0
        }
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            1
        }
    }
}
