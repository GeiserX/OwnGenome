//! The chromosome map of step 38: for each sequenced child with a
//! sequenced parent, the stretches of each chromosome the child got from
//! that parent, and, with a sequenced grandparent, which of the parent's
//! two copies each stretch is.
//!
//! A site tells something when the parent is heterozygous and the
//! genotypes give the allele the parent passed on. Stretches are smoothed:
//! a run of sites pointing to the other copy is kept only with
//! [`MIN_SITES`] sites over [`MIN_BP`] bases; a shorter run (a genotype
//! error, a gene conversion) is absorbed and counted as discordant. A
//! stretch with no site for [`MAX_GAP`] bases is split there. Autosomes,
//! and chrX outside the pseudoautosomal regions from the mother only.
//!
//! Writes, per child, `<child>_segments.tsv`, `<child>_dnapainter.csv`
//! (DNA Painter's import format, GRCh38 positions) and
//! `<child>_segments.svg`.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::path::Path;

use indexmap::IndexMap;
use og_pedigree::{Die, Members, Par, read_members};
use pycompat::fmt::fixed;

use super::args::{Args, Kind, Spec};
use super::vcf::{self, Alleles, Calls, resolve, set};
use super::{Result, print_out, write_file};

pub const SPEC: Spec = Spec {
    prog: "og family segments",
    usage: "usage: og family segments [-h] --members MEMBERS --vcf VCF [--par-bed PAR_BED] --out-dir OUT_DIR\n",
    options: &[
        ("--members", Kind::One),
        ("--vcf", Kind::One),
        ("--par-bed", Kind::One),
        ("--out-dir", Kind::One),
    ],
    required: &["--members", "--vcf", "--out-dir"],
};

/// A run of sites pointing to the other grandparent's copy needs this
/// much support to be a segment of its own.
pub const MIN_SITES: usize = 20;
pub const MIN_BP: i64 = 500_000;
/// No site for this many bases: the stretch is split.
pub const MAX_GAP: i64 = 5_000_000;
/// chr1 to chr22, then chrX: the map's chromosomes, in its order.
pub const N_CHROMS: usize = 23;
const X: usize = 22;

/// The index of a map chromosome.
pub fn chrom_index(c: &[u8]) -> Option<usize> {
    let n = c.strip_prefix(b"chr")?;
    if n == b"X" {
        return Some(X);
    }
    if n.is_empty() || n.len() > 2 || n[0] == b'0' || !n.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let v: usize = std::str::from_utf8(n).ok()?.parse().ok()?;
    (1..=22).contains(&v).then(|| v - 1)
}

fn chrom_name(i: usize) -> String {
    if i == X {
        "chrX".to_owned()
    } else {
        format!("chr{}", i + 1)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Side {
    Paternal,
    Maternal,
}

impl Side {
    fn as_str(self) -> &'static str {
        match self {
            Side::Paternal => "paternal",
            Side::Maternal => "maternal",
        }
    }
    fn other(self) -> Side {
        match self {
            Side::Paternal => Side::Maternal,
            Side::Maternal => Side::Paternal,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Basis {
    Grandparents,
    Parent,
    OtherCopy,
}

impl Basis {
    fn as_str(self) -> &'static str {
        match self {
            Basis::Grandparents => "parent and grandparents",
            Basis::Parent => "parent",
            Basis::OtherCopy => "other copy",
        }
    }
}

/// One child and one sequenced parent: the parent's informative sites.
struct Link {
    side: Side,
    parent: String,
    pcol: usize,
    gf: Option<String>,
    gm: Option<String>,
    gfcol: Option<usize>,
    gmcol: Option<usize>,
    /// Per chromosome, the positions where the allele passed on is known.
    sites: Vec<Vec<i64>>,
    /// Per chromosome, (position, 0 from the grandfather | 1 from the grandmother).
    hap: Vec<Vec<(i64, u8)>>,
}

impl Link {
    fn grandparents(&self) -> bool {
        self.gfcol.is_some() || self.gmcol.is_some()
    }

    fn ancestor(&self, label: u8) -> String {
        let name = if label == 0 { &self.gf } else { &self.gm };
        name.clone().unwrap_or_else(|| {
            format!(
                "{}'s {}",
                self.parent,
                if label == 0 { "father" } else { "mother" }
            )
        })
    }
}

struct ChildLinks {
    name: String,
    male: bool,
    col: usize,
    /// The child's father's column, when sequenced (chrX of a daughter).
    fcol: Option<usize>,
    /// The other parent's column per link side, when sequenced.
    other: [Option<usize>; 2],
    links: Vec<Link>,
}

/// The chromosome map being built from the joint VCF, record by record.
pub struct Collector {
    children: Vec<ChildLinks>,
    par: Par,
    sorted: vcf::Sorted,
    /// `##contig` lengths of the VCF's header.
    pub lengths: HashMap<String, i64>,
}

/// Whether a sequenced child has a sequenced parent: without one there
/// is no map and the VCF need not be read.
pub fn has_links(members: &Members) -> bool {
    let seq = |p: &Option<String>| {
        p.as_deref()
            .and_then(|p| members.get(p))
            .is_some_and(|m| m.sequenced)
    };
    members
        .values()
        .any(|m| m.sequenced && (seq(&m.father) || seq(&m.mother)))
}

/// The sequenced members, sorted.
pub fn sequenced(members: &Members) -> Vec<String> {
    let mut seq: Vec<String> = members
        .iter()
        .filter(|(_, m)| m.sequenced)
        .map(|(s, _)| s.clone())
        .collect();
    seq.sort();
    seq
}

/// `contig_lengths`: the `##contig=<ID=...,length=N>` lines of a header.
pub fn contig_length(line: &[u8]) -> Option<(String, i64)> {
    let rest = line.strip_prefix(b"##contig=<ID=")?;
    let end = rest
        .iter()
        .position(|&c| c == b',' || c == b'>')
        .unwrap_or(rest.len());
    if end == 0 || rest.get(end) != Some(&b',') {
        return None;
    }
    let id = String::from_utf8_lossy(&rest[..end]).into_owned();
    let mut tail = &rest[end + 1..];
    while let Some(k) = memchr::memmem::find(tail, b"length=") {
        let after = &tail[k + 7..];
        let n = after.iter().take_while(|c| c.is_ascii_digit()).count();
        if n > 0 {
            let v = std::str::from_utf8(&after[..n]).ok()?.parse().ok()?;
            return Some((id, v));
        }
        tail = &tail[k + 1..];
    }
    None
}

impl Collector {
    /// The links of every sequenced child and sequenced parent. `idx` maps
    /// a sample to its column.
    pub fn new(members: &Members, idx: &HashMap<String, usize>, par: Par) -> Collector {
        let col = |s: &Option<String>| -> Option<usize> {
            let s = s.as_deref()?;
            members.get(s).filter(|m| m.sequenced)?;
            idx.get(s).copied()
        };
        let mut children = Vec::new();
        for s in sequenced(members) {
            let m = &members[&s];
            let mut links = Vec::new();
            for (side, p) in [(Side::Paternal, &m.father), (Side::Maternal, &m.mother)] {
                let Some(pcol) = col(p) else { continue };
                let pm = &members[p.as_deref().unwrap_or_default()];
                links.push(Link {
                    side,
                    parent: p.clone().unwrap_or_default(),
                    pcol,
                    gf: pm.father.clone(),
                    gm: pm.mother.clone(),
                    gfcol: col(&pm.father),
                    gmcol: col(&pm.mother),
                    sites: vec![Vec::new(); N_CHROMS],
                    hap: vec![Vec::new(); N_CHROMS],
                });
            }
            if links.is_empty() {
                continue;
            }
            children.push(ChildLinks {
                male: m.sex.as_deref() == Some("male"),
                col: idx[&s],
                fcol: col(&m.father),
                other: [col(&m.mother), col(&m.father)],
                links,
                name: s,
            });
        }
        Collector {
            children,
            par,
            sorted: vcf::Sorted::default(),
            lengths: HashMap::new(),
        }
    }

    /// Whether any child has a link: without one the VCF need not be read.
    pub fn any(&self) -> bool {
        !self.children.is_empty()
    }

    /// A `##` header line.
    pub fn header(&mut self, line: &[u8]) {
        if let Some((id, n)) = contig_length(line) {
            self.lengths.insert(id, n);
        }
    }

    /// One record. `cols` are its columns, `pos` its parsed POS.
    pub fn record(
        &mut self,
        path: &str,
        pos: i64,
        calls: &mut Calls<'_>,
        cols: &[&[u8]],
        gi: Option<usize>,
    ) -> std::result::Result<(), Die> {
        if !self.any() {
            return Ok(());
        }
        self.sorted.next(path, cols[0], pos)?;
        let Some(ci) = chrom_index(cols[0]) else {
            return Ok(());
        };
        // Without a PAR file the pseudoautosomal sites cannot be told apart,
        // so chrX is left out.
        if ci == X && (self.par.is_empty() || self.par.contains("chrX", pos)) {
            return Ok(());
        }
        if gi.is_none() {
            return Ok(());
        }
        for ch in &mut self.children {
            let Some(cg) = calls.gt(ch.col)? else {
                continue;
            };
            for ln in &mut ch.links {
                let Some(pg) = calls.gt(ln.pcol)? else {
                    continue;
                };
                if pg.len() != 2 || pg[0] == pg[1] {
                    continue; // only a heterozygous parent tells which copy was passed on
                }
                let t = if ci == X {
                    if ln.side != Side::Maternal {
                        continue;
                    }
                    if ch.male {
                        if set(&cg).len() != 1 {
                            continue;
                        }
                        Some(cg[0])
                    } else {
                        let fg = match ch.fcol {
                            Some(c) => calls.gt(c)?,
                            None => None,
                        };
                        let Some(fg) = fg else { continue };
                        if set(&fg).len() != 1 || cg.len() != 2 {
                            continue;
                        }
                        passed_on(&cg, Some(&pg), Some(&[fg[0], fg[0]]), false)
                    }
                } else {
                    if cg.len() != 2 {
                        continue;
                    }
                    let k = usize::from(ln.side == Side::Maternal);
                    let other = match ch.other[k] {
                        Some(c) => calls.gt(c)?,
                        None => None,
                    };
                    passed_on(&cg, Some(&pg), other.as_deref(), ln.side == Side::Paternal)
                };
                let Some(t) = t.filter(|t| pg.contains(t)) else {
                    continue;
                };
                ln.sites[ci].push(pos);
                if !ln.grandparents() {
                    continue;
                }
                let gt_of = |c: Option<usize>,
                             calls: &mut Calls<'_>|
                 -> std::result::Result<Option<Alleles>, Die> {
                    match c {
                        Some(c) => calls.gt(c),
                        None => Ok(None),
                    }
                };
                // which of the parent's copies holds t: the parent's own origin here
                let label = if ci == X {
                    // the mother's X: one from her father (his only X), one from her mother
                    let gfg = gt_of(ln.gfcol, calls)?;
                    let gmg = gt_of(ln.gmcol, calls)?;
                    if gfg.as_ref().is_some_and(|g| set(g).len() != 1) {
                        continue;
                    }
                    let pair = gfg.as_ref().map(|g| [g[0], g[0]]);
                    let Some(po) = passed_on(
                        &pg,
                        gmg.as_deref(),
                        pair.as_ref().map(|p| p.as_slice()),
                        false,
                    ) else {
                        continue;
                    };
                    u8::from(t == po)
                } else {
                    let gf = gt_of(ln.gfcol, calls)?.map(|g| set(&g));
                    let gm = gt_of(ln.gmcol, calls)?.map(|g| set(&g));
                    let opts = resolve(&pg, gf.as_deref(), gm.as_deref());
                    let Some(from_gf) = unique(opts.iter().map(|o| o.0)) else {
                        continue;
                    };
                    u8::from(t != from_gf)
                };
                ln.hap[ci].push((pos, label));
            }
        }
        Ok(())
    }
}

/// The one value of an iterator, or `None` when it has none or several.
fn unique(mut it: impl Iterator<Item = u32>) -> Option<u32> {
    let first = it.next()?;
    it.all(|v| v == first).then_some(first)
}

/// The allele the parent passed on, when the genotypes give one.
fn passed_on(
    child: &[u32],
    parent: Option<&[u32]>,
    other: Option<&[u32]>,
    parent_is_father: bool,
) -> Option<u32> {
    let ps = parent.map(set);
    let os = other.map(set);
    let (f, m) = if parent_is_father { (ps, os) } else { (os, ps) };
    let opts = resolve(child, f.as_deref(), m.as_deref());
    unique(
        opts.iter()
            .map(|o| if parent_is_father { o.0 } else { o.1 }),
    )
}

/// Splits sorted positions where no site falls for [`MAX_GAP`] bases:
/// the index ranges of the blocks.
fn blocks(positions: impl Iterator<Item = i64>) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let (mut start, mut len) = (0, 0);
    let mut prev: Option<i64> = None;
    for (i, p) in positions.enumerate() {
        if let Some(q) = prev
            && p - q > MAX_GAP
        {
            out.push(start..i);
            start = i;
        }
        prev = Some(p);
        len = i + 1;
    }
    if len > start {
        out.push(start..len);
    }
    out
}

/// A run of one block: `[label, start, end, agree, discordant]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    pub label: u8,
    pub start: i64,
    pub end: i64,
    pub agree: usize,
    pub discordant: usize,
}

fn weak(r: &Run) -> bool {
    r.agree < MIN_SITES || r.end - r.start + 1 < MIN_BP
}

fn runs_of(sites: &[(i64, u8)]) -> Vec<Run> {
    let mut runs: Vec<Run> = Vec::new();
    for &(pos, label) in sites {
        match runs.last_mut() {
            Some(r) if r.label == label => {
                r.end = pos;
                r.agree += 1;
            }
            _ => runs.push(Run {
                label,
                start: pos,
                end: pos,
                agree: 1,
                discordant: 0,
            }),
        }
    }
    runs
}

/// Merges the weak run `i` with its neighbours `lo..=hi`: the merged run
/// takes the label of `runs[lo]` (or `runs[hi]` when `lo` is `i`), its
/// sites now point to the other copy, and the ones it had absorbed agree.
fn merge(lo: &Run, mid: &[&Run], hi: &Run, i_is: impl Fn(usize) -> bool, keep: u8) -> Run {
    let mut m = Run {
        label: keep,
        start: lo.start,
        end: hi.end,
        agree: 0,
        discordant: 0,
    };
    for (k, r) in mid.iter().enumerate() {
        if i_is(k) {
            m.agree += r.discordant;
            m.discordant += r.agree;
        } else {
            m.agree += r.agree;
            m.discordant += r.discordant;
        }
    }
    m
}

/// The tie-break of the weakest run: fewest sites, then shortest span,
/// then the first in the block (`order` grows along it).
type Key = (usize, i64, usize);

fn key(r: &Run, order: usize) -> Key {
    (r.agree, r.end - r.start, order)
}

/// The smoothed runs of one block of (position, label): repeatedly absorb
/// the weakest weak run into its neighbours. A min-heap with lazy
/// deletion over a linked list of runs, so O(R log R); the same merges, in
/// the same order, as taking the minimum of `(agree, end - start, index)`
/// over the list every time.
pub fn smooth(sites: &[(i64, u8)]) -> Vec<Run> {
    let runs = runs_of(sites);
    let n = runs.len();
    // Node i: its run, its list neighbours, its order (index of its first
    // original run), and whether it is still in the list.
    let mut node: Vec<Run> = runs;
    let mut prev: Vec<Option<usize>> = (0..n).map(|i| i.checked_sub(1)).collect();
    let mut next: Vec<Option<usize>> = (0..n).map(|i| (i + 1 < n).then_some(i + 1)).collect();
    let mut order: Vec<usize> = (0..n).collect();
    let mut alive: Vec<bool> = vec![true; n];
    let mut count = n;
    let mut heap: BinaryHeap<Reverse<(Key, usize)>> = BinaryHeap::new();
    for (i, r) in node.iter().enumerate() {
        if weak(r) {
            heap.push(Reverse((key(r, i), i)));
        }
    }
    while count > 1 {
        let Some(Reverse((_, i))) = heap.pop() else {
            break;
        };
        if !alive[i] {
            continue;
        }
        let lo = prev[i].unwrap_or(i);
        let hi = next[i].unwrap_or(i);
        let keep = if lo != i {
            node[lo].label
        } else {
            node[hi].label
        };
        let ids: Vec<usize> = [lo, i, hi].into_iter().fold(Vec::new(), |mut v, x| {
            if !v.contains(&x) {
                v.push(x);
            }
            v
        });
        let mid: Vec<&Run> = ids.iter().map(|&x| &node[x]).collect();
        let merged = merge(&node[lo], &mid, &node[hi], |k| ids[k] == i, keep);
        let id = node.len();
        node.push(merged);
        prev.push(prev[lo]);
        next.push(next[hi]);
        order.push(order[lo]);
        alive.push(true);
        if let Some(p) = prev[lo] {
            next[p] = Some(id);
        }
        if let Some(q) = next[hi] {
            prev[q] = Some(id);
        }
        for &x in &ids {
            alive[x] = false;
        }
        count -= ids.len() - 1;
        if weak(&node[id]) {
            heap.push(Reverse((key(&node[id], order[id]), id)));
        }
    }
    // Walk the list from its head.
    let mut out = Vec::new();
    let mut cur = (0..node.len()).find(|&x| alive[x] && prev[x].is_none());
    while let Some(x) = cur {
        if node[x].agree >= MIN_SITES {
            out.push(node[x]);
        }
        cur = next[x];
    }
    out
}

/// One row of `<child>_segments.tsv`.
#[derive(Clone, Debug)]
struct Row {
    side: Side,
    parent: String,
    chrom: usize,
    start: i64,
    end: i64,
    sites: usize,
    discordant: usize,
    ancestor: String,
    basis: Basis,
}

impl Row {
    fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// The map: per child (sorted), its rows in chromosome order, father's side first.
pub struct Map {
    rows: Vec<(String, Vec<Row>)>,
    /// `##contig` lengths of the VCF the map was read from.
    lengths: HashMap<String, i64>,
}

impl Collector {
    /// The rows of every child, from the sites collected.
    pub fn finish(self, members: &Members) -> Map {
        let lengths = self.lengths;
        let mut out = Vec::new();
        for ch in self.children {
            let mut rows = Vec::new();
            for ln in &ch.links {
                let total: usize = ln.hap.iter().map(Vec::len).sum();
                let hap = ln.grandparents() && total >= MIN_SITES;
                for c in 0..N_CHROMS {
                    if ln.sites[c].is_empty() {
                        continue;
                    }
                    if hap && !ln.hap[c].is_empty() {
                        let sites = &ln.hap[c];
                        for blk in blocks(sites.iter().map(|s| s.0)) {
                            for r in smooth(&sites[blk]) {
                                rows.push(Row {
                                    side: ln.side,
                                    parent: ln.parent.clone(),
                                    chrom: c,
                                    start: r.start,
                                    end: r.end,
                                    sites: r.agree,
                                    discordant: r.discordant,
                                    ancestor: ln.ancestor(r.label),
                                    basis: Basis::Grandparents,
                                });
                            }
                        }
                    } else {
                        let ps = &ln.sites[c];
                        for blk in blocks(ps.iter().copied()) {
                            if blk.len() >= MIN_SITES {
                                rows.push(Row {
                                    side: ln.side,
                                    parent: ln.parent.clone(),
                                    chrom: c,
                                    start: ps[blk.start],
                                    end: ps[blk.end - 1],
                                    sites: blk.len(),
                                    discordant: 0,
                                    ancestor: ".".to_owned(),
                                    basis: Basis::Parent,
                                });
                            }
                        }
                    }
                }
            }
            // One parent sequenced: the other copy over the same stretches.
            if let [ln] = ch.links.as_slice() {
                let other_side = ln.side.other();
                let m = &members[&ch.name];
                let other = match ln.side {
                    Side::Paternal => &m.mother,
                    Side::Maternal => &m.father,
                };
                for c in 0..N_CHROMS {
                    // the X a child got from the father is not seen from the mother's sites
                    if c == X || ln.sites[c].is_empty() {
                        continue;
                    }
                    let ps = &ln.sites[c];
                    for blk in blocks(ps.iter().copied()) {
                        if blk.len() >= MIN_SITES {
                            rows.push(Row {
                                side: other_side,
                                parent: other.clone().unwrap_or_else(|| ".".to_owned()),
                                chrom: c,
                                start: ps[blk.start],
                                end: ps[blk.end - 1],
                                sites: blk.len(),
                                discordant: 0,
                                ancestor: ".".to_owned(),
                                basis: Basis::OtherCopy,
                            });
                        }
                    }
                }
            }
            rows.sort_by_key(|r| (r.chrom, r.side != Side::Paternal, r.start));
            out.push((ch.name, rows));
        }
        Map { rows: out, lengths }
    }
}

const TSV: &str =
    "child\tside\tparent\tchrom\tstart\tend\tlength_bp\tsites\tdiscordant\tancestor\tbasis\n";
const DNAPAINTER: [&str; 12] = [
    "Chr",
    "Start",
    "End",
    "cM",
    "SNPs",
    "Match",
    "Confidence",
    "Group",
    "Side",
    "Notes",
    "Colour",
    "Ancestor",
];
const OTHER: &str = "#6b6b6b";

fn colours(side: Side) -> [&'static str; 2] {
    match side {
        Side::Paternal => ["#4f8fd6", "#9cc3ef"],
        Side::Maternal => ["#d65f8f", "#efa3c0"],
    }
}

fn side_only(side: Side) -> &'static str {
    match side {
        Side::Paternal => "#3b6fa8",
        Side::Maternal => "#a84870",
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// `x * scale` as `%.1f`, with `x` an int (Python: int times float).
fn f1(x: f64) -> String {
    fixed(x, 1)
}

fn svg(child: &str, rows: &[Row], lengths: &HashMap<String, i64>) -> String {
    let chroms: Vec<usize> = (0..N_CHROMS)
        .filter(|&c| lengths.contains_key(&chrom_name(c)) || rows.iter().any(|r| r.chrom == c))
        .collect();
    let longest = chroms
        .iter()
        .map(|&c| lengths.get(&chrom_name(c)).copied().unwrap_or(0))
        .chain(rows.iter().map(|r| r.end))
        .chain([1])
        .max()
        .unwrap_or(1);
    let (width, left, bar, row_h) = (900i64, 60i64, 9i64, 28i64);
    let scale = (width - left - 20) as f64 / longest as f64;
    let mut colour: IndexMap<(Side, String), &'static str> = IndexMap::new();
    for r in rows {
        if r.basis == Basis::OtherCopy {
            colour.insert((r.side, r.ancestor.clone()), OTHER);
        } else if r.ancestor == "." {
            colour.insert((r.side, ".".to_owned()), side_only(r.side));
        } else {
            let seen = colour
                .keys()
                .filter(|k| k.0 == r.side && k.1 != ".")
                .count();
            colour
                .entry((r.side, r.ancestor.clone()))
                .or_insert(colours(r.side)[seen % 2]);
        }
    }
    let height = 50 + row_h * chroms.len() as i64 + 20 * (colour.len() as i64 + 1);
    let mut l = vec![
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" font-family=\"sans-serif\" font-size=\"11\">"
        ),
        "<rect width=\"100%\" height=\"100%\" fill=\"#16181c\"/>".to_owned(),
        format!(
            "<text x=\"10\" y=\"20\" fill=\"#e6e6e6\" font-size=\"14\">{}: which copy of each chromosome came from whom (GRCh38)</text>",
            esc(child)
        ),
    ];
    for (i, &c) in chroms.iter().enumerate() {
        let y = 40 + i as i64 * row_h;
        let name = chrom_name(c);
        l.push(format!(
            "<text x=\"10\" y=\"{}\" fill=\"#e6e6e6\">{}</text>",
            y + bar + 3,
            esc(&name[3..])
        ));
        let n = lengths.get(&name).copied().unwrap_or_else(|| {
            rows.iter()
                .filter(|r| r.chrom == c)
                .map(|r| r.end)
                .chain([1])
                .max()
                .unwrap_or(1)
        });
        for k in 0..2 {
            l.push(format!(
                "<rect x=\"{left}\" y=\"{}\" width=\"{}\" height=\"{bar}\" fill=\"#2a2d33\"/>",
                y + k * (bar + 1),
                f1(n as f64 * scale)
            ));
        }
        for r in rows.iter().filter(|r| r.chrom == c) {
            let k = i64::from(r.side != Side::Paternal);
            let w = ((r.end - r.start) as f64 * scale).max(1.0);
            let who = if r.ancestor != "." {
                &r.ancestor
            } else {
                &r.parent
            };
            l.push(format!(
                "<rect class=\"seg\" x=\"{}\" y=\"{}\" width=\"{}\" height=\"{bar}\" fill=\"{}\"><title>{} {}:{}-{}, {}, {} sites</title></rect>",
                f1(left as f64 + r.start as f64 * scale),
                y + k * (bar + 1),
                f1(w),
                colour[&(r.side, r.ancestor.clone())],
                esc(r.side.as_str()),
                esc(&name),
                r.start,
                r.end,
                esc(who),
                r.sites
            ));
        }
    }
    let y = 50 + row_h * chroms.len() as i64;
    l.push(format!(
        "<text x=\"10\" y=\"{y}\" fill=\"#e6e6e6\">Upper bar: the copy from the father. Lower bar: from the mother. Grey: not mapped.</text>"
    ));
    let mut legend: Vec<(&(Side, String), &&str)> = colour.iter().collect();
    legend.sort_by(|a, b| (a.0.0.as_str(), &a.0.1).cmp(&(b.0.0.as_str(), &b.0.1)));
    for (j, ((side, anc), col)) in legend.into_iter().enumerate() {
        let label = if *col == OTHER {
            format!(
                "{} side, the copy not from the sequenced parent",
                side.as_str()
            )
        } else if anc == "." {
            format!("{} side, from the parent", side.as_str())
        } else {
            format!("{} side, from {anc}", side.as_str())
        };
        let j = j as i64;
        l.push(format!(
            "<rect x=\"10\" y=\"{}\" width=\"12\" height=\"10\" fill=\"{col}\"/><text x=\"28\" y=\"{}\" fill=\"#e6e6e6\">{}</text>",
            y + 10 + 20 * j,
            y + 19 + 20 * j,
            esc(&label)
        ));
    }
    l.push("</svg>".to_owned());
    let mut s = l.join("\n");
    s.push('\n');
    s
}

impl Map {
    /// Writes `<child>_segments.tsv`, `<child>_dnapainter.csv` and
    /// `<child>_segments.svg` for every child.
    pub fn write(&self, out_dir: &Path) -> std::result::Result<(), Die> {
        let lengths = &self.lengths;
        std::fs::create_dir_all(out_dir)
            .map_err(|e| Die::one(format!("cannot create {}: {e}", out_dir.display())))?;
        for (child, rows) in &self.rows {
            let mut t = String::from(TSV);
            for r in rows {
                t.push_str(&format!(
                    "{child}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    r.side.as_str(),
                    r.parent,
                    chrom_name(r.chrom),
                    r.start,
                    r.end,
                    r.length(),
                    r.sites,
                    r.discordant,
                    r.ancestor,
                    r.basis.as_str()
                ));
            }
            write_file(&out_dir.join(format!("{child}_segments.tsv")), t.as_bytes())?;
            let mut w = pycompat::csv::Writer::new(Vec::new(), ',');
            let csv_err =
                |e: std::io::Error| Die::one(format!("cannot write the DNA Painter CSV: {e}"));
            w.write_row(DNAPAINTER).map_err(csv_err)?;
            for r in rows {
                let who = if r.ancestor != "." {
                    r.ancestor.clone()
                } else if r.parent != "." {
                    r.parent.clone()
                } else {
                    format!("{} side", r.side.as_str())
                };
                let how = match r.basis {
                    Basis::Grandparents => "the parent's copy, told by a sequenced grandparent",
                    Basis::Parent => "from the sequenced parent",
                    Basis::OtherCopy => "the copy not from the sequenced parent",
                };
                let note = format!("GRCh38 positions; {how}; {} discordant sites", r.discordant);
                // DNA Painter reads Side and Ancestor from the first row of each Group, and Ancestor
                // takes only 1 (a known ancestor), so each grandparent (or parent) is a group of its own.
                let known = r.ancestor != "." || r.parent != ".";
                let name = chrom_name(r.chrom);
                w.write_row([
                    &name[3..],
                    &r.start.to_string(),
                    &r.end.to_string(),
                    "",
                    &r.sites.to_string(),
                    &format!("{child} from {who}"),
                    "",
                    &who,
                    r.side.as_str(),
                    &note,
                    "",
                    if known { "1" } else { "" },
                ])
                .map_err(csv_err)?;
            }
            write_file(
                &out_dir.join(format!("{child}_dnapainter.csv")),
                &w.into_inner(),
            )?;
            write_file(
                &out_dir.join(format!("{child}_segments.svg")),
                svg(child, rows, lengths).as_bytes(),
            )?;
        }
        Ok(())
    }

    /// The "Chromosome map" section of `family_report.md`.
    pub fn markdown(&self, members: &Members) -> Vec<String> {
        let mut l = vec![String::new(), "## Chromosome map".to_owned(), String::new()];
        if self.rows.is_empty() {
            l.push("No sequenced child has a sequenced parent, so there is no map.".to_owned());
            return l;
        }
        l.push(
            "For each child, the stretches of each chromosome that came from each parent, from the sites where \
             that parent is heterozygous and the genotypes give the allele passed on. When a parent's own parent \
             is sequenced too, each stretch also says which of the parent's two copies it is (the grandfather's or \
             the grandmother's), and a change between them is a crossover. `<child>_segments.tsv` lists them, \
             `<child>_dnapainter.csv` is the same list in DNA Painter's import format (GRCh38 positions), \
             `<child>_segments.svg` draws it."
                .to_owned(),
        );
        l.push(String::new());
        for (child, rows) in &self.rows {
            let m = &members[child];
            for (side, who) in [(Side::Paternal, "father"), (Side::Maternal, "mother")] {
                let mine: Vec<&Row> = rows.iter().filter(|r| r.side == side).collect();
                if mine.is_empty() {
                    continue;
                }
                let bp: i64 = mine.iter().map(|r| r.length()).sum();
                let mb = fixed(bp as f64 / 1e6, 1);
                let name = match side {
                    Side::Paternal => m.father.as_deref(),
                    Side::Maternal => m.mother.as_deref(),
                }
                .unwrap_or("not named");
                let n = mine.len();
                if mine.iter().any(|r| r.basis == Basis::OtherCopy) {
                    l.push(format!(
                        "- {child}: from the {who} ({name}, not sequenced): {n} stretches over {mb} Mb, the copy that did not \
                         come from the sequenced parent"
                    ));
                } else if mine.iter().any(|r| r.basis == Basis::Grandparents) {
                    let x = crossovers(rows, side);
                    let mut anc: Vec<&str> = mine.iter().map(|r| r.ancestor.as_str()).collect();
                    anc.sort_unstable();
                    anc.dedup();
                    let disc: usize = mine.iter().map(|r| r.discordant).sum();
                    l.push(format!(
                        "- {child}: from the {who} ({name}): {n} stretches over {mb} Mb, told apart as the copy from {}; \
                         {x} crossover{}, {disc} discordant sites absorbed",
                        anc.join(" or "),
                        if x == 1 { "" } else { "s" }
                    ));
                } else {
                    l.push(format!(
                        "- {child}: from the {who} ({name}): {n} stretches over {mb} Mb; no grandparent on this side is \
                         sequenced, so the two copies of the {who} are not told apart"
                    ));
                }
            }
        }
        l.push(String::new());
        l.push(
            "A stretch starts and ends at its first and last informative site: a crossover lies somewhere \
             between two stretches, not at a base. Short reads give no long-range phase, so without a sequenced \
             grandparent the map says only which parent a stretch came from."
                .to_owned(),
        );
        l
    }
}

/// Changes of ancestor between neighbouring segments of one chromosome.
fn crossovers(rows: &[Row], side: Side) -> usize {
    let mut n = 0;
    let mut prev: Option<&Row> = None;
    for r in rows {
        if r.side != side || r.basis != Basis::Grandparents {
            continue;
        }
        if let Some(p) = prev
            && p.chrom == r.chrom
            && p.ancestor != r.ancestor
        {
            n += 1;
        }
        prev = Some(r);
    }
    n
}

/// Reads the VCF's header: the `##` lines go to the collector, the
/// `#CHROM` line gives the sample columns. Returns the lines reader at
/// the first record.
pub fn read_header<R: std::io::BufRead>(
    lines: &mut og_vcf::Lines<R>,
    path: &str,
    mut on_meta: impl FnMut(&[u8]),
) -> std::result::Result<Option<Vec<String>>, Die> {
    let io = |e: std::io::Error| Die::io(path, &e);
    // Peek by reading; the caller handles the record that ends the header.
    loop {
        let Some(line) = lines.next_line().map_err(io)? else {
            return Ok(None);
        };
        if line.starts_with(b"##") {
            on_meta(line);
            continue;
        }
        if line.starts_with(b"#CHROM") {
            return Ok(Some(vcf::sample_names(line)));
        }
        return Err(Die::one(format!(
            "{path}: a record comes before the #CHROM line"
        )));
    }
}

/// Reads the whole `--vcf` for the map alone (`og family segments`, or
/// the second reader of `og family report --phased`).
pub fn collect(
    members: &Members,
    path: &str,
    par: Par,
) -> std::result::Result<(Collector, Option<Vec<String>>), Die> {
    let r = og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(|e| Die::io(path, &e))?;
    let mut lines = og_vcf::Lines::new(r);
    let mut lengths = HashMap::new();
    let names = read_header(&mut lines, path, |l| {
        if let Some((id, n)) = contig_length(l) {
            lengths.insert(id, n);
        }
    })?;
    let idx = names.as_deref().map(vcf::sample_index).unwrap_or_default();
    if !has_links(members) {
        let mut col = Collector::new(members, &idx, par);
        col.lengths = lengths;
        return Ok((col, names));
    }
    let Some(names) = names else {
        return Err(Die::one(format!("{path} has no #CHROM line")));
    };
    let lost: Vec<String> = sequenced(members)
        .into_iter()
        .filter(|s| !idx.contains_key(s))
        .collect();
    if !lost.is_empty() {
        return Err(Die::one(format!(
            "{path} has no column for {}",
            lost.join(", ")
        )));
    }
    let mut col = Collector::new(members, &idx, par);
    col.lengths = lengths;
    let mut n = lines.number();
    while let Some(line) = lines.next_line().map_err(|e| Die::io(path, &e))? {
        n += 1;
        let cols = vcf::columns(line);
        if cols.len() < 2 {
            return Err(Die::one(format!("{path} line {n}: not a VCF record")));
        }
        let pos = vcf::pos(cols[1], path, n)?;
        if chrom_index(cols[0]).is_none() {
            col.sorted.next(path, cols[0], pos)?;
            continue;
        }
        let Some(format) = cols.get(8) else {
            return Err(Die::one(format!("{path} line {n}: no FORMAT column")));
        };
        let gi = vcf::format_index(format, b"GT");
        let mut calls = Calls::new(&cols, gi.unwrap_or(0), path, n);
        col.record(path, pos, &mut calls, &cols, gi)?;
    }
    Ok((col, Some(names)))
}

pub fn run(a: &Args) -> Result<()> {
    let members = read_members(a.one("--members").unwrap_or_default())?;
    let vcf_path = a.one("--vcf").unwrap_or_default();
    let par = Par::read(a.one("--par-bed"))?;
    let (col, _) = collect(&members, vcf_path, par)?;
    let map = col.finish(&members);
    map.write(Path::new(a.one("--out-dir").unwrap_or_default()))?;
    let mut text = map.markdown(&members).join("\n");
    text.push('\n');
    print_out(&text);
    Ok(())
}

#[cfg(test)]
mod tests;
