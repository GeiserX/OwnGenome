//! `og y tree-positions` and `og y tree-place` (step 37): place a Y
//! chromosome on the current YFull YTree from the reads.
//!
//! ```text
//! og y tree-positions --tree-dir DIR OUT
//! og y tree-place --tree-dir DIR --pileup PILEUP --sample S --out-prefix P [--min-depth N]
//! ```
//!
//! Yleaf 4.1 places the sample on the YFull tree it ships (v14.01). This
//! tool walks the newer YTree release `versions.env` pins (`YTREE_VERSION`)
//! itself, so step 37 reports both placements, each with its tree.
//!
//! The tree folder (`setup.sh --yleaf-data` installs it as
//! `<yleaf data>/ytree`): `current_tree.json` (every branch with its id and
//! its defining SNPs by name, `A/B` being one SNP with two names), the two
//! GRCh38 SNP tables of Yleaf 4.x (`positions_ftdna_hg38.txt`,
//! `positions_yfull_hg38.txt`: name, position, ancestral and derived base;
//! a name is looked up in the first table, then the second) and
//! `version.txt` (`tree<TAB>...` and `positions<TAB>...` lines).
//!
//! `tree-positions` writes `chrY<TAB>pos` for every tree SNP the tables
//! place outside the pseudoautosomal regions, sorted, for
//! `samtools mpileup -l`. `tree-place` reads that pileup (made without
//! `-f`) and writes `P_y_tree.txt` (one result row) and `P_y_tree_path.tsv`
//! (one row per branch from the root to the terminal one).
//!
//! The rules. A SNP is called when at least `--min-depth` reads (10) cover
//! it and 90% of them carry one of its two bases: derived (D) or ancestral
//! (A); anything else is a no-call (N). A branch is positive with at least
//! one D and more D than A, negative with an A and no more D than A. The
//! terminal branch is the deepest positive branch with no blocking branch
//! above it (ties: more D on the path, then the smaller id). A negative
//! branch blocks the branches below it, except one with a single A: the
//! branches below it on the path outvote that site when they carry D at
//! [`OUTVOTE`] or more distinct positions, and the site is reported
//! (`Outvoted_ancestral`). Two A on one branch still block, and so does a
//! second single-A branch on the same path.
//!
//! How it is built: the tree is read as a stream of JSON events straight
//! into an arena of branches in the order the Python walked them (parents
//! first, children in file order), with SNP names interned once. The
//! position table is one slot per name the tree uses. The outvote needs only
//! "two distinct positions or more", so each branch carries a three-state
//! summary instead of a set of positions.
//!
//! What changes on purpose: the pileup is keyed by contig and position.
//! Only `chrY` (or `Y`) lines are read; lines on any other contig are
//! counted in the step line, and a pileup with both `chrY` and `Y` is
//! refused (the Python keyed by position only, so two contigs would mix).
//! A branch whose `id`, `snps` or `children` is of a type the Python could
//! not use, or that repeats one of those keys, is refused with its place in
//! the tree (the Python crashed, or kept the last copy). Both outputs are
//! written whole or not at all.

use crate::y::{is_file, path_join};
use crate::yjson::{Event, Reader};
use crate::{Error, Outcome};
use og_filters::args::Spec;
use og_io::Utf8;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Write;

/// Yleaf's `-r` default.
pub const MIN_DEPTH: i64 = 10;
/// Yleaf's `-b` default (90%).
const MAJORITY: f64 = 0.9;
/// Distinct derived sites below a single-ancestral-SNP branch that outvote it.
pub const OUTVOTE: usize = 2;
/// GRCh38 chrY PAR1 is 10,001-2,781,479.
pub const PAR1_END: u64 = 2_781_479;
/// GRCh38 chrY PAR2 is 56,887,903-57,217,415.
pub const PAR2_START: u64 = 56_887_903;
pub const TABLES: [&str; 2] = ["positions_ftdna_hg38.txt", "positions_yfull_hg38.txt"];
const TREE: &str = "current_tree.json";
const NONE: u32 = u32::MAX;

/// A base of a SNP table, as an index into the pileup counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Base {
    A,
    C,
    G,
    T,
}

impl Base {
    fn from_text(s: &str) -> Option<Base> {
        match s {
            "A" => Some(Base::A),
            "C" => Some(Base::C),
            "G" => Some(Base::G),
            "T" => Some(Base::T),
            _ => None,
        }
    }

    fn from_char(c: char) -> Option<Base> {
        match c {
            'A' | 'a' => Some(Base::A),
            'C' | 'c' => Some(Base::C),
            'G' | 'g' => Some(Base::G),
            'T' | 't' => Some(Base::T),
            _ => None,
        }
    }
}

/// A SNP a table places: position, ancestral and derived base.
#[derive(Debug, Clone, Copy)]
struct Snp {
    pos: u32,
    anc: Base,
    der: Base,
}

/// What the reads say about one SNP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Call {
    Derived,
    Ancestral,
    /// A no-call (N).
    None,
}

/// One branch of the tree.
struct Node {
    /// `None` for a missing, empty or false id.
    id: Option<Box<str>>,
    /// The parent's index, [`NONE`] for the root.
    parent: u32,
    depth: u32,
    /// Its SNPs: a range of [`Tree::groups`].
    groups: (u32, u32),
}

/// The tree as the Python's `walk` listed it, with the names interned.
struct Tree {
    nodes: Vec<Node>,
    /// One SNP each: a range of `aliases`.
    groups: Vec<(u32, u32)>,
    aliases: Vec<u32>,
    /// Each SNP name once, with its id. Ids index the position table; a
    /// name is looked up from its id only for the few printed.
    names: FxHashMap<Box<str>, u32>,
}

impl Tree {
    fn group(&self, g: u32) -> &[u32] {
        let (a, b) = self.groups[g as usize];
        &self.aliases[a as usize..b as usize]
    }

    fn node_groups(&self, n: &Node) -> std::ops::Range<u32> {
        n.groups.0..n.groups.1
    }

    /// The names of `ids`, in one pass over the names.
    fn names_of(&self, ids: &[u32]) -> HashMap<u32, &str> {
        let want: std::collections::HashSet<u32> = ids.iter().copied().collect();
        self.names
            .iter()
            .filter(|(_, id)| want.contains(id))
            .map(|(n, &id)| (id, n.as_ref()))
            .collect()
    }
}

/// Builds the arena from the JSON events.
struct Builder {
    tree: Tree,
}

impl Builder {
    fn intern(&mut self, name: &str) -> u32 {
        if let Some(&i) = self.tree.names.get(name) {
            return i;
        }
        let i = self.tree.names.len() as u32;
        self.tree.names.insert(name.into(), i);
        i
    }

    /// Python's `snp_groups`: split on `,`, each part on `/`, names stripped,
    /// empty names and groups dropped.
    fn set_snps(&mut self, node: usize, snps: &str) {
        let start = self.tree.groups.len() as u32;
        for g in snps.split(',') {
            let a = self.tree.aliases.len() as u32;
            for name in g.split('/') {
                let name = crate::strip(name);
                if !name.is_empty() {
                    let id = self.intern(name);
                    self.tree.aliases.push(id);
                }
            }
            let b = self.tree.aliases.len() as u32;
            if b > a {
                self.tree.groups.push((a, b));
            }
        }
        self.tree.nodes[node].groups = (start, self.tree.groups.len() as u32);
    }
}

enum Frame {
    /// Inside a branch's object; which of its keys were seen (id, snps, children).
    Node(usize, [bool; 3]),
    /// Inside a `children` list.
    Children(usize),
}

fn read_tree(text: &str) -> Result<Tree, String> {
    let mut b = Builder {
        tree: Tree {
            nodes: Vec::new(),
            groups: Vec::new(),
            aliases: Vec::new(),
            names: FxHashMap::default(),
        },
    };
    let mut r = Reader::new(text);
    let mut stack: Vec<Frame> = Vec::new();
    let new_node = |b: &mut Builder, parent: u32, depth: u32| {
        b.tree.nodes.push(Node {
            id: None,
            parent,
            depth,
            groups: (0, 0),
        });
        b.tree.nodes.len() - 1
    };
    match r.next()? {
        Event::StartObject => {
            let i = new_node(&mut b, NONE, 0);
            stack.push(Frame::Node(i, [false; 3]));
        }
        _ => return Err("the tree is not a JSON object".into()),
    }
    while let Some(top) = stack.last_mut() {
        match top {
            Frame::Node(i, seen) => {
                let i = *i;
                let key = match r.next()? {
                    Event::EndObject => {
                        stack.pop();
                        continue;
                    }
                    Event::Key(k) => k,
                    e => unreachable!("an object holds keys, not {e:?}"),
                };
                let slot = match key.as_ref() {
                    "id" => 0,
                    "snps" => 1,
                    "children" => 2,
                    _ => {
                        skip_value(&mut r)?;
                        continue;
                    }
                };
                if std::mem::replace(&mut seen[slot], true) {
                    return Err(format!(
                        "branch {} has two '{key}' keys",
                        where_(&b.tree, i)
                    ));
                }
                let v = r.next()?;
                match (slot, v) {
                    (0 | 1, Event::Str(s)) if !s.is_empty() => {
                        if slot == 0 {
                            b.tree.nodes[i].id = Some(s.as_ref().into());
                        } else {
                            b.set_snps(i, &s);
                        }
                    }
                    (2, Event::StartArray) => {
                        stack.push(Frame::Children(i));
                    }
                    (2, Event::StartObject) if r.next()? == Event::EndObject => {}
                    (_, v) if v.falsy() => {}
                    _ => {
                        return Err(format!(
                            "branch {}: '{key}' is not {}",
                            where_(&b.tree, i),
                            if slot == 2 { "a list" } else { "text" }
                        ));
                    }
                }
            }
            Frame::Children(parent) => {
                let parent = *parent;
                match r.next()? {
                    Event::EndArray => {
                        stack.pop();
                    }
                    Event::StartObject => {
                        let depth = b.tree.nodes[parent].depth + 1;
                        let i = new_node(&mut b, parent as u32, depth);
                        stack.push(Frame::Node(i, [false; 3]));
                    }
                    _ => {
                        return Err(format!(
                            "a child of branch {} is not an object",
                            where_(&b.tree, parent)
                        ));
                    }
                }
            }
        }
    }
    if r.next()? != Event::End {
        return Err("extra data after the tree".into());
    }
    Ok(b.tree)
}

/// A branch named for a message: its id, or its place in the walk.
fn where_(tree: &Tree, i: usize) -> String {
    match &tree.nodes[i].id {
        Some(id) => id.to_string(),
        None => format!("#{i} (no id yet)"),
    }
}

/// Reads past one value whose start is the next event.
fn skip_value(r: &mut Reader<'_>) -> Result<(), String> {
    let mut depth = 0usize;
    loop {
        match r.next()? {
            Event::StartObject | Event::StartArray => depth += 1,
            Event::EndObject | Event::EndArray => depth -= 1,
            Event::Key(_) => continue,
            _ => {}
        }
        if depth == 0 {
            return Ok(());
        }
    }
}

/// The tree folder read: the tree, the table slot of each name, the labels.
struct Loaded {
    tree: Tree,
    table: Vec<Option<Snp>>,
    labels: HashMap<String, String>,
}

fn load(dir: &str) -> Result<Loaded, Error> {
    let tree_path = path_join(dir, TREE);
    for need in [
        tree_path.clone(),
        path_join(dir, TABLES[0]),
        path_join(dir, TABLES[1]),
    ] {
        if !is_file(&need) {
            return Err(Error(format!(
                "{need} not found: install the current Y tree with scripts/setup.sh --yleaf-data <genome_dir>"
            )));
        }
    }
    let bytes = std::fs::read(&tree_path).map_err(|e| Error::io(&e, &tree_path))?;
    let text =
        std::str::from_utf8(&bytes).map_err(|e| Error(format!("{tree_path}: not UTF-8: {e}")))?;
    let tree = read_tree(text).map_err(|e| Error(format!("{tree_path}: {e}")))?;
    drop(bytes);
    // One slot per name the tree uses; the first valid row of the first
    // table that has the name fills it.
    let mut table: Vec<Option<Snp>> = vec![None; tree.names.len()];
    for t in TABLES {
        let path = path_join(dir, t);
        for line in crate::lines(&path, Utf8::Strict)? {
            let line = line?;
            let c: Vec<&str> = line
                .strip_suffix('\n')
                .unwrap_or(&line)
                .split('\t')
                .collect();
            if c.len() < 7 {
                continue;
            }
            let Some(&id) = tree.names.get(c[1]) else {
                continue;
            };
            if table[id as usize].is_some() || !crate::digits(c[3]) {
                continue;
            }
            let (Some(anc), Some(der)) = (Base::from_text(c[5]), Base::from_text(c[6])) else {
                continue;
            };
            if anc == der {
                continue;
            }
            // Too many digits for u64 is past PAR2 too.
            let pos = c[3].parse::<u64>().unwrap_or(u64::MAX);
            if pos <= PAR1_END || pos >= PAR2_START {
                continue;
            }
            table[id as usize] = Some(Snp {
                pos: pos as u32,
                anc,
                der,
            });
        }
    }
    let mut labels = HashMap::new();
    let vf = path_join(dir, "version.txt");
    if is_file(&vf) {
        let text = std::fs::read(&vf).map_err(|e| Error::io(&e, &vf))?;
        let text = String::from_utf8(text).map_err(|e| Error(format!("{vf}: not UTF-8: {e}")))?;
        for line in splitlines(&text) {
            let (k, v) = line.split_once('\t').unwrap_or((line, ""));
            labels.insert(crate::strip(k).to_owned(), crate::strip(v).to_owned());
        }
    }
    Ok(Loaded {
        tree,
        table,
        labels,
    })
}

/// Python's `str.splitlines()`.
fn splitlines(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut start = 0;
    let mut it = s.char_indices().peekable();
    while let Some((i, c)) = it.next() {
        let brk = matches!(
            c,
            '\n' | '\r'
                | '\u{b}'
                | '\u{c}'
                | '\u{1c}'
                | '\u{1d}'
                | '\u{1e}'
                | '\u{85}'
                | '\u{2028}'
                | '\u{2029}'
        );
        if brk {
            out.push(&s[start..i]);
            let mut next = i + c.len_utf8();
            if c == '\r' && it.peek().map(|&(_, n)| n) == Some('\n') {
                it.next();
                next += 1;
            }
            start = next;
        }
    }
    if start < s.len() {
        out.push(&s[start..]);
    }
    out
}

/// `tree-positions`: the sorted distinct positions; the step line.
fn positions(dir: &str, out: &str) -> Outcome {
    let l = load(dir).map_err(|e| (e, 1))?;
    let mut pos: Vec<u32> = l.table.iter().flatten().map(|s| s.pos).collect();
    pos.sort_unstable();
    pos.dedup();
    if pos.is_empty() {
        return Err((
            Error(format!("no tree SNP is on the position tables in {dir}")),
            1,
        ));
    }
    let mut body = String::with_capacity(pos.len() * 14);
    for p in &pos {
        let _ = writeln!(body, "chrY\t{p}");
    }
    crate::write_atomic(out, body.as_bytes()).map_err(|e| (e, 1))?;
    Ok(format!(
        "{} Y positions of the current tree ({} branches)\n",
        pos.len(),
        l.tree.nodes.len()
    ))
}

/// Counts of A, C, G, T in a samtools mpileup read-base column: `^` skips
/// the mapping-quality character after it, `+N`/`-N` skip the N inserted or
/// deleted bases, `.` and `,` are the reference base, the rest of the
/// characters that are not bases (`$`, `*`, `<`, `>`) count for nothing.
fn base_counts(bases: &str, reference: &str) -> [u32; 4] {
    let mut counts = [0u32; 4];
    let reference = match reference {
        "A" | "a" => Some(Base::A),
        "C" | "c" => Some(Base::C),
        "G" | "g" => Some(Base::G),
        "T" | "t" => Some(Base::T),
        _ => None,
    };
    let mut it = bases.chars();
    while let Some(ch) = it.next() {
        let b = match ch {
            '^' => {
                it.next();
                continue;
            }
            '+' | '-' => {
                let rest = it.as_str();
                let n = rest.bytes().take_while(u8::is_ascii_digit).count();
                if n > 0 {
                    // Skip the digits, then that many bases.
                    let skip: usize = rest[..n].parse().unwrap_or(usize::MAX);
                    it = rest[n..].chars();
                    for _ in 0..skip {
                        if it.next().is_none() {
                            break;
                        }
                    }
                }
                continue;
            }
            '.' | ',' => reference,
            c => Base::from_char(c),
        };
        if let Some(b) = b {
            counts[b as usize] += 1;
        }
    }
    counts
}

/// The pileup: counts by position on the Y contig; lines on other contigs.
struct Pileup {
    counts: HashMap<u32, [u32; 4]>,
    other_contigs: u64,
}

fn read_pileup(path: &str) -> Result<Pileup, Error> {
    let mut counts = HashMap::new();
    let mut other = 0u64;
    let mut seen = [false; 2];
    for line in crate::lines(path, Utf8::Strict)? {
        let line = line?;
        let c: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(&line)
            .split('\t')
            .collect();
        if c.len() < 5 || !crate::digits(c[1]) {
            continue;
        }
        match c[0] {
            "chrY" => seen[0] = true,
            "Y" => seen[1] = true,
            _ => {
                other += 1;
                continue;
            }
        }
        if seen == [true, true] {
            return Err(Error(format!(
                "{path} holds both chrY and Y lines: make the pileup on one Y contig"
            )));
        }
        // A position past u32 is on no table.
        if let Ok(pos) = c[1].parse::<u32>() {
            counts.insert(pos, base_counts(c[4], c[2]));
        }
    }
    Ok(Pileup {
        counts,
        other_contigs: other,
    })
}

fn call(snp: Snp, pileup: &Pileup, min_depth: i128) -> Call {
    let Some(c) = pileup.counts.get(&snp.pos) else {
        return Call::None;
    };
    let depth: u64 = c.iter().map(|&x| u64::from(x)).sum();
    if i128::from(depth) < min_depth {
        return Call::None;
    }
    // As the Python compares: the count against 0.9 times the depth, in floats.
    let limit = MAJORITY * depth as f64;
    if f64::from(c[snp.der as usize]) >= limit {
        Call::Derived
    } else if f64::from(c[snp.anc as usize]) >= limit {
        Call::Ancestral
    } else {
        Call::None
    }
}

/// The distinct positions of derived SNPs seen so far, as far as the
/// outvote needs them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Sites {
    Zero,
    One(u32),
    TwoOrMore,
}

impl Sites {
    fn add(self, pos: u32) -> Sites {
        match self {
            Sites::Zero => Sites::One(pos),
            Sites::One(p) if p == pos => self,
            _ => Sites::TwoOrMore,
        }
    }

    fn union(self, other: Sites) -> Sites {
        match other {
            Sites::Zero => self,
            Sites::One(p) => self.add(p),
            Sites::TwoOrMore => Sites::TwoOrMore,
        }
    }

    fn count(self) -> usize {
        match self {
            Sites::Zero => 0,
            Sites::One(_) => 1,
            Sites::TwoOrMore => 2,
        }
    }
}

/// A branch's calls.
#[derive(Debug, Clone, Copy)]
struct Stats {
    snps: u32,
    on_table: u32,
    d: u32,
    a: u32,
    n: u32,
    der_sites: Sites,
}

/// Python's `place_snp`: the first alias a table places, or the first alias.
fn place_snp(group: &[u32], table: &[Option<Snp>]) -> (u32, Option<Snp>) {
    group
        .iter()
        .find_map(|&n| table[n as usize].map(|s| (n, Some(s))))
        .unwrap_or((group[0], None))
}

fn node_stats(l: &Loaded, node: &Node, pileup: &Pileup, min_depth: i128) -> Stats {
    let mut s = Stats {
        snps: 0,
        on_table: 0,
        d: 0,
        a: 0,
        n: 0,
        der_sites: Sites::Zero,
    };
    for g in l.tree.node_groups(node) {
        s.snps += 1;
        let (_, Some(snp)) = place_snp(l.tree.group(g), &l.table) else {
            continue;
        };
        s.on_table += 1;
        match call(snp, pileup, min_depth) {
            Call::Derived => {
                s.d += 1;
                s.der_sites = s.der_sites.add(snp.pos);
            }
            Call::Ancestral => s.a += 1,
            Call::None => s.n += 1,
        }
    }
    s
}

/// The name ids of a branch's ancestral SNPs, in order.
fn ancestral_ids(l: &Loaded, node: &Node, pileup: &Pileup, min_depth: i128) -> Vec<u32> {
    l.tree
        .node_groups(node)
        .filter_map(|g| match place_snp(l.tree.group(g), &l.table) {
            (name, Some(snp)) if call(snp, pileup, min_depth) == Call::Ancestral => Some(name),
            _ => None,
        })
        .collect()
}

/// `tree-place`: the two files and the step line.
fn place(dir: &str, pileup_path: &str, sample: &str, prefix: &str, min_depth: i128) -> Outcome {
    let fail = |e| (e, 1);
    let l = load(dir).map_err(fail)?;
    let pileup = read_pileup(pileup_path).map_err(fail)?;
    let nodes = &l.tree.nodes;
    let stats: Vec<Stats> = nodes
        .iter()
        .map(|n| node_stats(&l, n, &pileup, min_depth))
        .collect();
    let positive = |s: &Stats| s.d >= 1 && s.d > s.a;
    let negative = |s: &Stats| s.a >= 1 && s.d <= s.a;
    let single = |s: &Stats| negative(s) && s.a == 1;
    // Down the walk (parents come first): the D on the path, the
    // single-ancestral-SNP branches on it, whether a blocking branch is on
    // it, and the derived sites below its single-ancestral-SNP branch (None
    // without one).
    let n = nodes.len();
    let mut path_d = vec![0u64; n];
    let mut singles = vec![0u8; n];
    let mut blocked = vec![false; n];
    let mut soft: Vec<Option<Sites>> = vec![None; n];
    let mut best: Option<usize> = None;
    for (i, node) in nodes.iter().enumerate() {
        let s = &stats[i];
        let p = (node.parent != NONE).then_some(node.parent as usize);
        path_d[i] = u64::from(s.d) + p.map_or(0, |p| path_d[p]);
        singles[i] = u8::from(single(s)).saturating_add(p.map_or(0, |p| singles[p]));
        blocked[i] =
            (negative(s) && !single(s)) || singles[i] >= 2 || p.is_some_and(|p| blocked[p]);
        soft[i] = if single(s) {
            Some(Sites::Zero)
        } else {
            p.and_then(|p| soft[p]).map(|x| x.union(s.der_sites))
        };
        if positive(s) && !blocked[i] && soft[i].is_none_or(|x| x.count() >= OUTVOTE) {
            let better = match best {
                None => true,
                Some(b) => {
                    let (nb, nn) = (&nodes[b], node);
                    let key = (nn.depth, path_d[i]);
                    let bkey = (nb.depth, path_d[b]);
                    key > bkey || (key == bkey && id_or_empty(nn) < id_or_empty(nb))
                }
            };
            if better {
                best = Some(i);
            }
        }
    }
    // The path from the root to the best branch, branches without an id left out.
    let mut path = Vec::new();
    let mut at = best.map_or(NONE, |b| b as u32);
    while at != NONE {
        let node = &nodes[at as usize];
        if node.id.is_some() {
            path.push(at as usize);
        }
        at = node.parent;
    }
    path.reverse();
    let tree_label = l.labels.get("tree").map_or("current tree", String::as_str);
    let pos_label = l
        .labels
        .get("positions")
        .map_or("position tables", String::as_str);
    let called: u64 = stats.iter().map(|s| u64::from(s.d + s.a)).sum();
    let term = path.last().copied();
    let id = |i: usize| nodes[i].id.as_deref().unwrap_or_default();
    let hg = term.map_or("NA", id);
    let mut body = String::from(
        "Level\tBranch\tSNPs\tOn_tables\tDerived\tAncestral\tNo_call\tAncestral_SNPs\n",
    );
    let mut anc_term = String::from(".");
    let mut outvoted: Vec<&str> = Vec::new();
    let (mut pd, mut pa, mut pn) = (0u64, 0u64, 0u64);
    let anc_ids: Vec<Vec<u32>> = path
        .iter()
        .map(|&i| ancestral_ids(&l, &nodes[i], &pileup, min_depth))
        .collect();
    let names = l.tree.names_of(&anc_ids.concat());
    for (lvl, (&i, ids)) in path.iter().zip(&anc_ids).enumerate() {
        let s = &stats[i];
        let anc: Vec<&str> = ids.iter().map(|id| names[id]).collect();
        let joined = anc.join(",");
        let _ = writeln!(
            body,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            lvl + 1,
            id(i),
            s.snps,
            s.on_table,
            s.d,
            s.a,
            s.n,
            if joined.is_empty() { "." } else { &joined }
        );
        (pd, pa, pn) = (
            pd + u64::from(s.d),
            pa + u64::from(s.a),
            pn + u64::from(s.n),
        );
        if single(s) {
            outvoted.push(anc[0]);
        }
        if Some(i) == term && !joined.is_empty() {
            anc_term = joined;
        }
    }
    let outvoted = if outvoted.is_empty() {
        ".".to_owned()
    } else {
        outvoted.join(",")
    };
    let level = path.len();
    let row = format!(
        "Sample_name\tHg\tLevel\tPath_derived\tPath_ancestral\tPath_no_call\tAncestral_in_terminal\t\
         SNPs_called\tMin_depth\tTree\tPositions\tOutvoted_ancestral\n\
         {sample}\t{hg}\t{level}\t{pd}\t{pa}\t{pn}\t{anc_term}\t{called}\t{min_depth}\t\
         {tree_label}\t{pos_label}\t{outvoted}\n"
    );
    crate::write_atomic(&format!("{prefix}_y_tree_path.tsv"), body.as_bytes()).map_err(fail)?;
    crate::write_atomic(&format!("{prefix}_y_tree.txt"), row.as_bytes()).map_err(fail)?;
    let mut line = if term.is_none() {
        format!(
            "Current tree ({tree_label}): no branch has a derived call ({called} SNPs called at depth >= {min_depth})"
        )
    } else {
        let mut extra = String::new();
        if anc_term != "." {
            let _ = write!(extra, "; ancestral inside its block: {anc_term}");
        }
        if outvoted != "." {
            let _ = write!(
                extra,
                "; ancestral on the path, outvoted by the SNPs below it: {outvoted}"
            );
        }
        format!(
            "Current tree ({tree_label}): {hg} (level {level}, path SNPs {pd} derived, {pa} ancestral, {pn} no call{extra})"
        )
    };
    if pileup.other_contigs > 0 {
        let _ = write!(
            line,
            " [{} pileup {} on another contig than chrY, ignored]",
            pileup.other_contigs,
            if pileup.other_contigs == 1 {
                "line"
            } else {
                "lines"
            }
        );
    }
    line.push('\n');
    Ok(line)
}

fn id_or_empty(n: &Node) -> &str {
    n.id.as_deref().unwrap_or_default()
}

const POSITIONS: Spec = Spec {
    prog: "og y tree-positions",
    usage: "usage: og y tree-positions [-h] --tree-dir TREE_DIR out\n",
    options: &["--tree-dir"],
    positionals: &["out"],
    required: &["--tree-dir"],
};

const PLACE: Spec = Spec {
    prog: "og y tree-place",
    usage: "usage: og y tree-place [-h] --tree-dir TREE_DIR --pileup PILEUP --sample SAMPLE --out-prefix OUT_PREFIX [--min-depth MIN_DEPTH]\n",
    options: &[
        "--tree-dir",
        "--pileup",
        "--sample",
        "--out-prefix",
        "--min-depth",
    ],
    positionals: &[],
    required: &["--tree-dir", "--pileup", "--sample", "--out-prefix"],
};

pub(crate) fn positions_cli(rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    crate::run(&POSITIONS, rest, out, err, |a| {
        Ok(positions(crate::opt(a, "--tree-dir"), &a.pos[0]))
    })
}

pub(crate) fn place_cli(rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    crate::run(&PLACE, rest, out, err, |a| {
        let min_depth = a.typed(
            "--min-depth",
            i128::from(MIN_DEPTH),
            "int",
            og_filters::args::int,
        )?;
        if min_depth < 1 {
            return Ok(Err((Error("--min-depth must be 1 or more".into()), 1)));
        }
        Ok(place(
            crate::opt(a, "--tree-dir"),
            crate::opt(a, "--pileup"),
            crate::opt(a, "--sample"),
            crate::opt(a, "--out-prefix"),
            min_depth,
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    #[test]
    fn pileup_marks() {
        // Read starts with their quality, ends, indels; '.'/',' as the reference.
        assert_eq!(base_counts("^I..,$+2ACt-1Ga*", "t"), [1, 0, 0, 4]);
        assert_eq!(base_counts("+x.", "A"), [1, 0, 0, 0]);
        assert_eq!(base_counts("..", "N"), [0, 0, 0, 0]);
        assert_eq!(base_counts("^", "A"), [0, 0, 0, 0]);
        assert_eq!(base_counts("+12AC", "A"), [0, 0, 0, 0]);
        assert_eq!(base_counts("<>acgt", "A"), [1, 1, 1, 1]);
    }

    #[test]
    fn splitlines_as_python() {
        assert_eq!(splitlines("a\r\nb\rc\u{85}d\n"), ["a", "b", "c", "d"]);
        assert_eq!(splitlines("a\n\nb"), ["a", "", "b"]);
        assert!(splitlines("").is_empty());
    }

    #[test]
    fn the_walk_order_is_preorder_in_file_order() {
        let t = read_tree(
            r#"{"children": [{"id": "A", "children": [{"id": "C"}]}, {"id": "B"}], "snps": "", "id": ""}"#,
        )
        .unwrap();
        let ids: Vec<_> = t
            .nodes
            .iter()
            .map(|n| n.id.as_deref().unwrap_or("-"))
            .collect();
        assert_eq!(ids, ["-", "A", "C", "B"]);
        assert_eq!(t.nodes[2].parent, 1);
        assert_eq!(t.nodes[2].depth, 2);
    }

    #[test]
    fn snps_split_into_alias_groups() {
        let t =
            read_tree(r#"{"children": [{"id": "A", "snps": " a1, a2/ a2alt ,, / ,b"}]}"#).unwrap();
        let groups: Vec<Vec<&str>> = t
            .node_groups(&t.nodes[1])
            .map(|g| {
                let ids = t.group(g);
                let names = t.names_of(ids);
                ids.iter().map(|n| names[n]).collect()
            })
            .collect();
        assert_eq!(groups, [vec!["a1"], vec!["a2", "a2alt"], vec!["b"]]);
    }

    #[test]
    fn falsy_and_refused_branch_values() {
        assert!(read_tree(r#"{"id": null, "snps": 0, "children": {}}"#).is_ok());
        assert!(
            read_tree(r#"{"id": false, "snps": "", "children": null, "x": [1, {"y": []}]}"#)
                .is_ok()
        );
        for bad in [
            r#"{"id": 5}"#,
            r#"{"snps": ["a"]}"#,
            r#"{"children": "ab"}"#,
            r#"{"children": [1]}"#,
            r#"{"children": {"a": 1}}"#,
            r#"{"id": "a", "id": "b"}"#,
            r#"[]"#,
            r#"{"children": []} {}"#,
        ] {
            assert!(read_tree(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn outvote_sites() {
        assert_eq!(Sites::Zero.add(5).add(5), Sites::One(5));
        assert_eq!(Sites::One(5).union(Sites::One(6)).count(), 2);
        assert_eq!(Sites::One(5).union(Sites::One(5)).count(), 1);
        assert_eq!(Sites::Zero.union(Sites::TwoOrMore).count(), 2);
    }

    /// The tree of tests/test_y_tree_place.py, small: root -> A (a1) -> X (x1)
    /// -> Y (y1, y1b at one site) -> Z (z1).
    fn folder(d: &Dir) {
        d.put(
            "current_tree.json",
            r#"{"id": "", "children": [{"id": "A", "snps": "a1", "children": [
               {"id": "X", "snps": "x1", "children": [{"id": "Y", "snps": "y1, y1b", "children": [
               {"id": "Z", "snps": "z1"}]}]}]}]}"#,
        );
        d.put(
            "positions_ftdna_hg38.txt",
            "chry\ta1\tR\t3000100\tC->T\tC\tT\nchry\tx1\tR\t3000200\tA->C\tA\tC\n\
             chry\ty1\tR\t3000300\tT->A\tT\tA\nchry\ty1b\tR\t3000300\tT->A\tT\tA\n\
             chry\tz1\tR\t3000400\tG->C\tG\tC\n",
        );
        d.put("positions_yfull_hg38.txt", "");
    }

    fn pile(d: &Dir, rows: &[(u32, &str)]) -> String {
        let body: String = rows
            .iter()
            .map(|(p, b)| format!("chrY\t{p}\tN\t12\t{}\t{}\n", b.repeat(12), "I".repeat(12)))
            .collect();
        d.put("p.txt", &body)
    }

    fn hg(d: &Dir) -> String {
        let t = d.read("o_y_tree.txt");
        t.lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .nth(1)
            .unwrap()
            .to_owned()
    }

    /// Red first: with OUTVOTE 1, one derived site below x1 would outvote it.
    #[test]
    fn red_first_two_names_at_one_site_do_not_outvote() {
        let d = Dir::new("ytree-outvote");
        folder(&d);
        let dir = d.root().display().to_string();
        let prefix = d.path("o");
        let p = pile(&d, &[(3000100, "T"), (3000200, "A"), (3000300, "A")]);
        place(&dir, &p, "S", &prefix, 10).unwrap();
        assert_eq!(
            hg(&d),
            "A",
            "x1 ancestral, one site below (two names): blocked"
        );
        let p = pile(
            &d,
            &[
                (3000100, "T"),
                (3000200, "A"),
                (3000300, "A"),
                (3000400, "C"),
            ],
        );
        place(&dir, &p, "S", &prefix, 10).unwrap();
        assert_eq!(hg(&d), "Z", "two distinct sites below outvote x1");
        assert!(d.read("o_y_tree.txt").ends_with("\tx1\n"));
    }

    #[test]
    fn pileup_contigs() {
        let d = Dir::new("ytree-contigs");
        folder(&d);
        let dir = d.root().display().to_string();
        let p = d.put(
            "p.txt",
            "chrY\t3000100\tN\t12\tTTTTTTTTTTTT\tI\nchr1\t3000100\tN\t12\tCCCCCCCCCCCC\tI\n",
        );
        let line = place(&dir, &p, "S", &d.path("o"), 10).unwrap();
        assert_eq!(hg(&d), "A", "the chr1 line does not overwrite the chrY one");
        assert!(
            line.ends_with("[1 pileup line on another contig than chrY, ignored]\n"),
            "{line}"
        );
        let p = d.put(
            "p.txt",
            "chrY\t3000100\tN\t1\tT\tI\nY\t3000100\tN\t1\tT\tI\n",
        );
        assert!(place(&dir, &p, "S", &d.path("o"), 10).is_err());
    }
}
