//! Verify SURVIVOR clusters against the exact VCFs passed to the merge.
//! See docs/22-survivor-merge.md for the event and output rules.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};

use og_filters::args::{self, ArgError, Spec};

use super::info_value;

const USAGE: &str = "usage: og sv consensus-check --callers LIST IN.vcf[.gz]|- OUT.vcf|-
";
const SPEC: Spec = Spec {
    prog: "og sv consensus-check",
    usage: USAGE,
    options: &["--callers"],
    positionals: &["IN", "OUT"],
    required: &["--callers"],
};

const MIN_OVERLAP: f64 = 0.5;
const INS_MAX_DIST: u64 = 1000;
const INS_MIN_RATIO: f64 = 0.5;
const BND_MAX_DIST: u64 = 1000;
const MAX_GROUP_COMBINATIONS: usize = 1_000_000;

pub const HEADER: &str = "##og_sv_consensus_check=<Description=\"Each record has two or more caller records verified against the input VCFs: same type and strands, reciprocal overlap >= 0.5; insertions within 1000 bp and size ratio >= 0.5; intercontig breakends agree on both contigs and are within 1000 bp at both ends. Unresolved members do not vote.\">";

const EMPTY: &str = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
const SURVIVOR_FORMAT: &str = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO";

#[derive(Debug, Clone)]
pub(super) struct Member {
    pub caller: usize,
    pub kind: String,
    pub strands: String,
    pub chrom: String,
    pub chrom2: String,
    pub start: i64,
    pub end: i64,
    pub len: i64,
    pub id: String,
    pub reference: String,
    pub alt: String,
    pub qual: String,
    pub gt: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Apart {
    Type,
    Strands,
    Overlap,
    Distance,
    Size,
}

fn type_class(kind: &str, intercontig: bool) -> Option<&'static str> {
    match kind.split(':').next()? {
        "DEL" => Some("DEL"),
        "DUP" => Some("DUP"),
        "INV" => Some("INV"),
        "TRA" => Some("TRA"),
        "INS" | "ALU" | "LINE1" | "SVA" => Some("INS"),
        "BND" if intercontig => Some("TRA"),
        "BND" => Some("INV"),
        "CNV" => Some("CNV"),
        _ => None,
    }
}

fn numeric_qual(qual: &str) -> Option<f64> {
    let value: f64 = qual.parse().ok()?;
    value.is_finite().then_some(value)
}

pub(super) fn apart(a: &Member, b: &Member) -> Option<Apart> {
    let kind = type_class(&a.kind, a.chrom != a.chrom2);
    if kind != type_class(&b.kind, b.chrom != b.chrom2) {
        return Some(Apart::Type);
    }
    if a.strands != b.strands {
        return Some(Apart::Strands);
    }
    if a.chrom != b.chrom || a.chrom2 != b.chrom2 {
        return Some(Apart::Distance);
    }
    if kind == Some("INS") {
        if a.start.abs_diff(b.start) > INS_MAX_DIST {
            return Some(Apart::Distance);
        }
        let long = a.len.max(b.len);
        return (long == 0 || (a.len.min(b.len) as f64) < INS_MIN_RATIO * long as f64)
            .then_some(Apart::Size);
    }
    if a.chrom != a.chrom2 {
        return (a.start.abs_diff(b.start) > BND_MAX_DIST || a.end.abs_diff(b.end) > BND_MAX_DIST)
            .then_some(Apart::Distance);
    }
    let (a0, a1) = (a.start.min(a.end), a.start.max(a.end));
    let (b0, b1) = (b.start.min(b.end), b.start.max(b.end));
    let overlap = a1.min(b1) - a0.max(b0);
    let longer = (a1 - a0).max(b1 - b0);
    (overlap <= 0 || (overlap as f64) < MIN_OVERLAP * longer as f64).then_some(Apart::Overlap)
}

/// Search caller by caller. First choice wins a tie, so groups of equal
/// size are ordered by caller column and then by the member's CO index.
fn largest_group(
    remaining: &[Vec<usize>],
    members: &[&Member],
    group: &mut Vec<usize>,
    best: &mut Vec<usize>,
) {
    if group.len() > best.len() {
        best.clone_from(group);
    }
    if group.len() + remaining.len() <= best.len() {
        return;
    }
    if let Some((choices, rest)) = remaining.split_first() {
        for &i in choices {
            if group
                .iter()
                .all(|&j| apart(members[j], members[i]).is_none())
            {
                group.push(i);
                largest_group(rest, members, group, best);
                group.pop();
                if group.len() + remaining.len() <= best.len() {
                    return;
                }
            }
        }
        largest_group(rest, members, group, best);
    }
}

/// Visit caller columns and their CO entries in order. Add each call to
/// the first compatible group with room for its caller.
fn greedy_groups(candidates: Vec<Vec<usize>>, members: &[&Member]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for choices in candidates {
        for i in choices {
            let group = groups.iter_mut().find(|group| {
                group.iter().all(|&j| {
                    members[j].caller != members[i].caller
                        && apart(members[j], members[i]).is_none()
                })
            });
            if let Some(group) = group {
                group.push(i);
            } else {
                groups.push(vec![i]);
            }
        }
    }
    groups.retain(|group| group.len() >= 2);
    groups
}

/// Assign the largest compatible group first within the search budget.
/// Remove its members and repeat while two callers can still agree.
fn groups(members: &[&Member], caller_count: usize, counts: &mut Counts) -> Vec<Vec<usize>> {
    let mut candidates = vec![Vec::new(); caller_count];
    for (i, member) in members.iter().enumerate() {
        candidates[member.caller].push(i);
    }
    // Each column also has a branch for skipping that caller.
    let combinations = candidates.iter().fold(1_usize, |total, choices| {
        total.saturating_mul(choices.len() + 1)
    });
    if combinations > MAX_GROUP_COMBINATIONS {
        counts.greedy_clusters += 1;
        return greedy_groups(candidates, members);
    }
    let mut groups = Vec::new();
    loop {
        let mut best = Vec::new();
        largest_group(&candidates, members, &mut Vec::new(), &mut best);
        if best.len() < 2 {
            return groups;
        }
        for &i in &best {
            candidates[members[i].caller].retain(|&j| j != i);
        }
        groups.push(best);
    }
}

/// The pinned merge writes these numeric coordinates and a type for every
/// member. It prints the cluster's second contig even for a member with a
/// different mate. Use only the caller, first contig, numeric coordinates
/// and type to locate the input row. Every collision is ambiguous.
type Key = (String, i64, i64, String);

struct Caller {
    name: String,
    members: HashMap<Key, Option<Member>>,
}

fn position(value: &str) -> Option<i64> {
    let pos: i32 = value.parse().ok()?;
    (pos > 0).then_some(i64::from(pos))
}

fn allele_len(s: &str) -> Option<usize> {
    (!s.is_empty() && s.bytes().all(|b| b.is_ascii_alphabetic())).then_some(s.len())
}

/// Both brackets must agree. The local bases identify which side is joined.
fn breakend(alt: &str) -> Option<(&str, i64, &'static str)> {
    let i = alt.find(['[', ']'])?;
    let bracket = alt.as_bytes()[i] as char;
    let (mate, tail) = alt[i + 1..].split_once(bracket)?;
    let (chrom, pos) = mate.rsplit_once(':')?;
    let left = &alt[..i];
    let strands = match (allele_len(left), allele_len(tail), bracket) {
        (Some(_), None, ']') if tail.is_empty() => "++",
        (Some(_), None, '[') if tail.is_empty() => "+-",
        (None, Some(_), ']') if left.is_empty() => "-+",
        (None, Some(_), '[') if left.is_empty() => "--",
        _ => return None,
    };
    (!chrom.is_empty()).then_some((chrom, position(pos)?, strands))
}

fn source_strands(kind: &str, info: &str, mate: Option<(&str, i64, &str)>) -> Option<String> {
    if let Some((_, _, strands)) = mate {
        return Some(strands.to_owned());
    }
    let strands = info_value(info, "STRANDS");
    if !matches!(strands, "" | ".") {
        let mut pairs = strands.split(',').map(|s| s.split(':').next().unwrap_or(s));
        let pair = pairs.next()?;
        return (matches!(pair, "++" | "+-" | "-+" | "--") && pairs.all(|p| p == pair))
            .then(|| pair.to_owned());
    }
    let pair = match info_value(info, "CT") {
        "3to3" | "NtoN" => "++",
        "3to5" => "+-",
        "5to3" => "-+",
        "5to5" => "--",
        "" | "." => match kind {
            "DEL" | "INS" => "+-",
            "DUP" => "-+",
            "INV" | "CNV" => "++",
            _ => return None,
        },
        _ => return None,
    };
    Some(pair.to_owned())
}

/// SURVIVOR uses atoi, which stops at the first nondecimal character.
fn survivor_integer(value: &str) -> i64 {
    let value = value.trim_start();
    let unsigned = value.strip_prefix(['-', '+']).unwrap_or(value);
    let digits = unsigned.bytes().take_while(u8::is_ascii_digit).count();
    let end = value.len() - unsigned.len() + digits;
    value[..end].parse::<i32>().map_or(0, i64::from)
}

/// Mirror the pinned merge parser so its CO coordinates can join to the
/// prepared calls. The first recognized length in INFO wins.
fn survivor_length(info: &str) -> Option<i64> {
    for (i, field) in info.split(';').enumerate() {
        if let Some(value) = field.strip_prefix("SVLEN=") {
            // The first-field case in SURVIVOR skips two extra bytes.
            let value = if i == 0 {
                value.get(2..).unwrap_or_default()
            } else {
                value
            };
            return Some(survivor_integer(value).abs());
        }
        if i > 0
            && let Some(value) = field.strip_prefix("AVGLEN=")
        {
            let value = value.split(',').next().unwrap_or_default();
            return Some(i64::from(value.parse::<f64>().unwrap_or_default() as i32).abs());
        }
        // The INSLEN guard in 1.0.7 reads from the '=' and never succeeds.
    }
    None
}

fn source_member(line: &str, caller: usize) -> Option<(Key, Option<Member>)> {
    let fields: Vec<_> = line.split('\t').collect();
    let [chrom, start, id, reference, alt, qual, _, info, rest @ ..] = fields.as_slice() else {
        return None;
    };
    let (chrom, id, reference, alt, qual, info) = (*chrom, *id, *reference, *alt, *qual, *info);
    let start = position(start)?;
    let mate = breakend(alt);
    let chrom2 = mate.map_or_else(
        || match info_value(info, "CHR2") {
            "" | "." => chrom,
            other => other,
        },
        |(other, _, _)| other,
    );
    let symbolic = alt.strip_prefix('<').and_then(|s| s.strip_suffix('>'));
    let sequence_len = allele_len(reference)
        .zip(allele_len(alt))
        .map(|(a, b)| b as i64 - a as i64);
    let kind = match info_value(info, "SVTYPE") {
        "" | "." => {
            if let Some(allele) = symbolic {
                allele
            } else if mate.is_some() {
                "BND"
            } else if sequence_len? > 0 {
                "INS"
            } else {
                "DEL"
            }
        }
        kind => kind,
    };
    let kind = type_class(kind, chrom != chrom2)?;
    // SURVIVOR's END parser recognizes the semicolon-prefixed key. With
    // no recognized END, CO uses POS + length, including for insertions
    // (the merge adds length to the end and subtracts it on output).
    let stated_len = survivor_length(info);
    let ins_len = stated_len.or(sequence_len).filter(|&n| n > 0);
    let end = mate
        .map(|(_, end, _)| end)
        .or_else(|| {
            info.rsplit_once(";END=")
                .map(|(_, value)| survivor_integer(value))
                .filter(|&end| end != -1)
        })
        .or_else(|| Some(start + stated_len.or(sequence_len.map(i64::abs))?))?;
    let end = if kind == "INS" { end.max(1) } else { end };
    let key = (chrom.to_owned(), start, end, kind.to_owned());

    let len = match kind {
        "INS" => ins_len,
        _ if chrom != chrom2 => Some(0),
        _ => Some((end - start).abs()),
    };
    // Keep the join key even when a source call is unusable, so another
    // source row with the same key cannot silently inherit its vote.
    let member = len
        .zip(source_strands(kind, info, mate))
        .and_then(|(len, strands)| {
            if (alt.contains(['[', ']']) && mate.is_none())
                || (chrom != chrom2 && kind != "TRA")
                || (symbolic.is_some_and(|alt_type| {
                    type_class(alt_type, chrom != chrom2).is_some_and(|alt_kind| alt_kind != kind)
                }))
            {
                return None;
            }
            let gt = rest.first().zip(rest.get(1)).and_then(|(format, sample)| {
                let i = format.split(':').position(|name| name == "GT")?;
                sample.split(':').nth(i)
            });
            Some(Member {
                caller,
                kind: kind.to_owned(),
                strands,
                chrom: chrom.to_owned(),
                chrom2: chrom2.to_owned(),
                start,
                end,
                len,
                id: id.to_owned(),
                reference: reference.to_owned(),
                alt: alt.to_owned(),
                qual: qual.to_owned(),
                gt: gt.unwrap_or("./.").to_owned(),
            })
        });
    Some((key, member))
}

fn invalid(path: &str, message: &str) -> Error {
    Error::Read(
        path.to_owned(),
        io::Error::new(io::ErrorKind::InvalidData, message),
    )
}

impl Caller {
    fn read(reader: impl BufRead, path: &str, caller: usize) -> Result<Self, Error> {
        let mut name = None;
        let mut members = HashMap::new();
        for line in reader.lines() {
            let line = line.map_err(|e| Error::Read(path.to_owned(), e))?;
            if line.starts_with("#CHROM\t") {
                let columns: Vec<_> = line.split('\t').skip(9).collect();
                if let [sample] = columns.as_slice()
                    && !sample.is_empty()
                    && name.is_none()
                {
                    name = Some((*sample).to_owned());
                } else {
                    return Err(invalid(path, "expected one caller column"));
                }
            } else if !line.starts_with('#')
                && let Some((key, member)) = source_member(&line, caller)
            {
                members
                    .entry(key)
                    .and_modify(|value| *value = None)
                    .or_insert(member);
            }
        }
        Ok(Self {
            name: name.ok_or_else(|| invalid(path, "missing caller header"))?,
            members,
        })
    }
}

fn absent(v: &str) -> bool {
    matches!(v, "" | "." | "NAN" | "NaN" | "NA")
}

/// CO is CHROM_POS-CHROM2_POS2. Both contig names can contain '_' and '-' .
fn join_key(entry: &str, chrom: &str, kind: &str) -> Option<Key> {
    let rest = entry.strip_prefix(chrom)?.strip_prefix('_')?;
    let (start, mate) = rest.split_once('-')?;
    let (other, end) = mate.rsplit_once('_')?;
    if other.is_empty() {
        return None;
    }
    Some((
        chrom.to_owned(),
        position(start)?,
        position(end)?,
        kind.to_owned(),
    ))
}

struct Record<'a> {
    members: Vec<&'a Member>,
    unresolved: usize,
}

impl<'a> Record<'a> {
    fn read(line: &str, callers: &'a [Caller]) -> Option<Self> {
        let fields: Vec<_> = line.split('\t').collect();
        let [chrom, _, _, _, _, _, _, info, format, columns @ ..] = fields.as_slice() else {
            return None;
        };
        let bits = info_value(info, "SUPP_VEC").as_bytes();
        let supp = bits.iter().filter(|&&b| b == b'1').count();
        if supp < 2
            || bits.len() != callers.len()
            || columns.len() != callers.len()
            || bits.iter().any(|&b| !matches!(b, b'0' | b'1'))
            || info_value(info, "SUPP").parse() != Ok(supp)
        {
            return None;
        }
        let at = |name: &str| format.split(':').position(|key| key == name);
        let (co, ty) = (at("CO")?, at("TY")?);
        let ln = at("LN");
        let mut members = Vec::new();
        let mut unresolved = 0;
        for ((caller, column), &bit) in callers.iter().zip(columns).zip(bits) {
            let vals: Vec<_> = column.split(':').collect();
            let value = |i| vals.get(i).copied().unwrap_or(".");
            let coords = value(co);
            if bit == b'0' {
                if !absent(coords) {
                    return None;
                }
                continue;
            }
            if absent(coords) {
                unresolved += 1;
                continue;
            }
            let coords: Vec<_> = coords.split(',').collect();
            let kinds: Vec<_> = value(ty).split(',').collect();
            if coords.len() != kinds.len() {
                unresolved += coords.len();
                continue;
            }
            let insertion_vote_ok = ln.and_then(|i| position(value(i))).is_some();
            let mut seen = HashSet::new();
            for (entry, kind) in coords.into_iter().zip(kinds) {
                let member = join_key(entry, chrom, kind)
                    .filter(|_| kind != "INS" || insertion_vote_ok)
                    .and_then(|key| caller.members.get(&key))
                    .and_then(Option::as_ref)
                    .filter(|member| seen.insert((member.start, member.end, member.kind.as_str())));
                if let Some(member) = member {
                    members.push(member);
                } else {
                    unresolved += 1;
                }
            }
        }
        Some(Self {
            members,
            unresolved,
        })
    }
}

fn format_member(member: &Member) -> String {
    let id = member.id.replace(':', "_");
    let alt = if allele_len(&member.alt).is_some() {
        &member.alt
    } else {
        "NA"
    };
    format!(
        "{}:NA:{}:.,.:{}:{}:{}:{}:{}:{}:{}_{}-{}_{}",
        member.gt,
        member.len,
        member.strands,
        member.qual,
        member.kind,
        id,
        member.reference,
        alt,
        member.chrom,
        member.start,
        member.chrom2,
        member.end
    )
}

fn build_record(kept: &[&Member], caller_count: usize) -> String {
    let mut ordered = kept.to_vec();
    ordered.sort_by_key(|m| (m.start, m.end));
    let mid = ordered[(kept.len() - 1) / 2];
    // Keep an original bracketed allele, with its own POS and REF, when
    // some caller used a symbolic ALT for the same intercontig breakend.
    let mid = if mid.chrom != mid.chrom2 && breakend(&mid.alt).is_none() {
        ordered
            .iter()
            .copied()
            .find(|m| breakend(&m.alt).is_some())
            .unwrap_or(mid)
    } else {
        mid
    };
    let named = |id: &str| !matches!(id, "" | ".");
    let id = if named(&mid.id) {
        &mid.id
    } else {
        kept.iter()
            .map(|member| &member.id)
            .find(|id| named(id))
            .map_or(".", String::as_str)
    };
    let max_qual = kept
        .iter()
        .filter_map(|m| numeric_qual(&m.qual).map(|n| (n, m.qual.as_str())))
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .map_or(".", |(_, qual)| qual);
    let mut svlen =
        (kept.iter().map(|m| m.len).sum::<i64>() as f64 / kept.len() as f64).round() as i64;
    if mid.kind == "DEL" {
        svlen = -svlen;
    }

    let mut bits = vec!['0'; caller_count];
    let mut columns = vec![EMPTY.to_owned(); caller_count];
    for member in kept {
        bits[member.caller] = '1';
        columns[member.caller] = format_member(member);
    }
    let bits: String = bits.into_iter().collect();
    let (pos, end) = (mid.start, mid.end);
    let ci = |coord: fn(&Member) -> i64, point: i64| {
        let (mut min, mut max) = (point, point);
        for member in kept {
            min = min.min(coord(member));
            max = max.max(coord(member));
        }
        format!("{},{}", min - point, max - point)
    };
    format!(
        "{}\t{pos}\t{}\t{}\t{}\t{max_qual}\tPASS\tSUPP={};SUPP_VEC={bits};SVLEN={svlen};SVTYPE={};SVMETHOD=SURVIVOR1.0.7;CHR2={};END={end};CIPOS={};CIEND={};STRANDS={}\t{SURVIVOR_FORMAT}\t{}",
        mid.chrom,
        id,
        mid.reference,
        mid.alt,
        kept.len(),
        mid.kind,
        mid.chrom2,
        ci(|m| m.start, pos),
        ci(|m| m.end, end),
        mid.strands,
        columns.join("\t")
    )
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(super) struct Counts {
    pub records: usize,
    pub kept: usize,
    pub split: usize,
    pub from_split: usize,
    pub dropped: usize,
    pub unread: usize,
    pub unresolved: usize,
    pub greedy_clusters: usize,
}

impl Counts {
    pub fn written(&self) -> usize {
        self.kept + self.from_split
    }
}

impl fmt::Display for Counts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} records, {} written; {} one event as merged, {} rechecked: split into {} records of 2+ callers, {} dropped with fewer than 2 callers left; {} malformed records dropped, {} unresolved member votes dropped",
            self.records,
            self.written(),
            self.kept,
            self.split,
            self.from_split,
            self.dropped,
            self.unread,
            self.unresolved
        )?;
        if self.greedy_clusters > 0 {
            write!(
                f,
                "; clusters using greedy grouping: {}",
                self.greedy_clusters
            )?;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub(super) enum Error {
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

fn check_line(line: &str, callers: &[Caller], counts: &mut Counts) -> Vec<String> {
    counts.records += 1;
    let Some(record) = Record::read(line, callers) else {
        counts.unread += 1;
        return Vec::new();
    };
    counts.unresolved += record.unresolved;
    let groups = groups(&record.members, callers.len(), counts);
    if let [group] = groups.as_slice()
        && record.unresolved == 0
        && group.len() == record.members.len()
    {
        counts.kept += 1;
        return vec![line.to_owned()];
    }
    counts.split += 1;
    if groups.is_empty() {
        counts.dropped += 1;
    }
    counts.from_split += groups.len();
    groups
        .iter()
        .map(|group| {
            let kept: Vec<_> = group.iter().map(|&i| record.members[i]).collect();
            build_record(&kept, callers.len())
        })
        .collect()
}

fn check(
    input: impl BufRead,
    path: &str,
    callers: &[Caller],
    out: &mut impl Write,
) -> Result<Counts, Error> {
    let w = |e| Error::Write("the output".to_owned(), e);
    let mut counts = Counts::default();
    let mut saw_header = false;
    for line in input.lines() {
        let line = line.map_err(|e| Error::Read(path.to_owned(), e))?;
        if line.is_empty() {
            continue;
        }
        if line.starts_with("#CHROM\t") {
            let names = line.split('\t').skip(9);
            if saw_header || names.ne(callers.iter().map(|caller| caller.name.as_str())) {
                return Err(invalid(
                    path,
                    "caller columns differ from the supplied VCF list",
                ));
            }
            saw_header = true;
            writeln!(out, "{HEADER}").map_err(w)?;
        }
        if line.starts_with('#') {
            writeln!(out, "{line}").map_err(w)?;
        } else {
            if !saw_header {
                return Err(invalid(path, "missing caller header"));
            }
            for rec in check_line(&line, callers, &mut counts) {
                writeln!(out, "{rec}").map_err(w)?;
            }
        }
    }
    if !saw_header {
        return Err(invalid(path, "missing caller header"));
    }
    out.flush().map_err(w)?;
    Ok(counts)
}

#[cfg(test)]
pub(super) fn check_text(
    input: &str,
    sources: &[&str],
    out: &mut impl Write,
) -> Result<Counts, Error> {
    let callers = sources
        .iter()
        .enumerate()
        .map(|(i, source)| Caller::read(source.as_bytes(), "caller.vcf", i))
        .collect::<Result<Vec<_>, _>>()?;
    check(input.as_bytes(), "merged.vcf", &callers, out)
}

fn run(list: &str, input: &str, out_path: &str, stdout: &mut dyn Write) -> Result<Counts, Error> {
    let paths = std::fs::read_to_string(list).map_err(|e| Error::Read(list.to_owned(), e))?;
    let mut callers = Vec::new();
    let mut names = HashSet::new();
    for path in paths.lines() {
        let reader = og_io::open(path).map_err(|e| Error::Read(path.to_owned(), e))?;
        let caller = Caller::read(reader, path, callers.len())?;
        if !names.insert(caller.name.clone()) {
            return Err(invalid(list, "duplicate caller names"));
        }
        callers.push(caller);
    }
    if callers.len() < 2 {
        return Err(invalid(list, "expected at least two caller VCFs"));
    }
    let reader = og_io::open(input).map_err(|e| Error::Read(input.to_owned(), e))?;
    if out_path == "-" {
        return check(reader, input, &callers, &mut BufWriter::new(stdout));
    }
    let part = format!("{out_path}.part");
    let file = File::create(&part).map_err(|e| Error::Write(part.clone(), e))?;
    let mut w = BufWriter::new(file);
    let counts = check(reader, input, &callers, &mut w).inspect_err(|_| {
        let _ = std::fs::remove_file(&part);
    })?;
    drop(w);
    std::fs::rename(&part, out_path).map_err(|e| Error::Write(out_path.to_owned(), e))?;
    Ok(counts)
}

/// Exit 0 on success, 1 on I/O errors, 2 on incorrect usage.
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
    match run(p.get("--callers").unwrap(), &p.pos[0], &p.pos[1], out) {
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
