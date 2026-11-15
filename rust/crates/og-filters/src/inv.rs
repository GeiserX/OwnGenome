//! `og sv inv-filter RAW CONVERTED [--pair-bp N] [--max-single-bp N]
//! [--short-bp N]`: keep a whole-span inversion only where Manta saw both
//! of its junctions (replaces `bin/manta_inv_filter.py`, which ran on the
//! Manta image's Python 2.7).
//!
//! Manta writes each inversion junction as a breakend pair flagged `INV3`
//! or `INV5`, and its `convertInversion.py` turns every same-chromosome
//! pair into one `SVTYPE=INV` record from POS to END, so a lone junction
//! whose breakends lie 150 Mb apart becomes a 150 Mb inversion. This writes
//! the converted file again, with:
//! - `INV_JUNCTIONS=2` on an INV record whose opposite junction is there
//!   (POS and END each within `--pair-bp`, or the same Manta `EVENT`);
//! - `INV_JUNCTIONS=1` on one without; one spanning `--max-single-bp` or
//!   more is replaced by its two breakends from Manta's own file;
//! - `INVDUP_CANDIDATE` on a single junction shorter than `--short-bp`.
//!
//! The records come out sorted by the `##contig` order, then POS; ties keep
//! the replacement breakends first, then the converted records in file
//! order (the Python's stable sort). Lines are read as Python 2.7 read
//! them: split on `\n` only, blank lines dropped.
//!
//! Fix against the Python: a record without the 8 VCF columns, or a POS
//! that is not a number, was a traceback (exit 1); now a message naming the
//! file and line, exit 2.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, BufRead, Write};

use crate::args::{self, ArgError, Spec};

const SPEC: Spec = Spec {
    prog: "og sv inv-filter",
    usage: "usage: og sv inv-filter [--pair-bp N] [--max-single-bp N] [--short-bp N] RAW_VCF CONVERTED_VCF > OUT.vcf\n",
    options: &["--pair-bp", "--max-single-bp", "--short-bp"],
    positionals: &["raw", "converted"],
    required: &[],
};

/// The three bounds, in bases.
#[derive(Debug, Clone, Copy)]
pub struct Bounds {
    pub pair_bp: i128,
    pub max_single_bp: i128,
    pub short_bp: i128,
}

impl Default for Bounds {
    fn default() -> Self {
        Bounds {
            pair_bp: 1000,
            max_single_bp: 1_000_000,
            short_bp: 1000,
        }
    }
}

/// Why the files could not be filtered.
#[derive(Debug)]
pub enum InvError {
    Io(String, io::Error),
    /// A record the filter cannot read: file, line (1-based) and why.
    Record(String, u64, String),
}

impl fmt::Display for InvError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InvError::Io(path, e) => write!(f, "{}", crate::os_error(e, path)),
            InvError::Record(path, line, why) => write!(f, "{path}:{line}: {why}"),
        }
    }
}

/// One VCF file as Python 2.7's line loop saw it: the `#` lines, and each
/// record as its line number and fields.
/// A record: its line number and its columns.
type Record = (u64, Vec<Vec<u8>>);

struct Vcf {
    header: Vec<Vec<u8>>,
    records: Vec<Record>,
}

fn read_vcf(path: &str) -> Result<Vcf, InvError> {
    let io_err = |e| InvError::Io(path.to_owned(), e);
    let mut input = og_io::open(path).map_err(io_err)?;
    let mut vcf = Vcf {
        header: Vec::new(),
        records: Vec::new(),
    };
    let mut line = Vec::new();
    let mut number = 0;
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line).map_err(io_err)? == 0 {
            break;
        }
        number += 1;
        while line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' {
            vcf.header.push(line.clone());
        } else {
            vcf.records.push((
                number,
                line.split(|&c| c == b'\t').map(<[u8]>::to_vec).collect(),
            ));
        }
    }
    Ok(vcf)
}

/// An INFO value: `key=value`, or a bare flag.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum Value {
    Str(Vec<u8>),
    Flag,
}

/// INFO as a dict, the last of a repeated key winning.
fn info_dict(info: &[u8]) -> HashMap<&[u8], Value> {
    let mut d = HashMap::new();
    if info.is_empty() || info == b"." {
        return d;
    }
    for item in info.split(|&c| c == b';') {
        match item.iter().position(|&c| c == b'=') {
            Some(k) => d.insert(&item[..k], Value::Str(item[k + 1..].to_vec())),
            None => d.insert(item, Value::Flag),
        };
    }
    d
}

fn add_info(fields: &mut [Vec<u8>], item: &[u8]) {
    let info = &mut fields[7];
    if info.is_empty() || info == b"." {
        *info = item.to_vec();
    } else {
        info.push(b';');
        info.extend_from_slice(item);
    }
}

/// `int()` of Python 2.7: whitespace, a sign and ASCII digits.
fn py2_int(s: &[u8]) -> Option<i128> {
    let t = s.trim_ascii();
    if t.contains(&b'_') {
        return None;
    }
    pycompat::num::parse_int(std::str::from_utf8(t).ok()?).ok()
}

fn bad_pos(pos: &[u8]) -> String {
    format!(
        "POS {} is not a whole number",
        crate::py_repr(&String::from_utf8_lossy(pos))
    )
}

/// `MantaBND:12:0:1:0:0:0:0` -> `:12:0:1:0:0:0`, the part
/// `convertInversion.py` keeps in the INV record's id after `MantaINV`.
fn pair_key(id: &[u8]) -> Option<&[u8]> {
    let rest = id.strip_prefix(b"MantaBND")?;
    let i = rest.iter().rposition(|&c| c == b':')?;
    (i > 0).then(|| &rest[..i])
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Three,
    Five,
}

/// One INV record of the converted file.
struct Inv {
    /// Its index among the records.
    n: usize,
    pos: i128,
    end: i128,
    side: Side,
    event: Option<Value>,
}

/// The records (by index) that have their opposite junction: an INV3 and
/// an INV5 record on one chromosome with POS and END each within
/// `pair_bp`, or with the same EVENT.
fn paired(invs: &[Inv], chroms: &[&[u8]], pair_bp: i128) -> Vec<bool> {
    let mut out = vec![false; chroms.len()];
    // INV5 records per chromosome, sorted by POS, for the window search;
    // and per (chromosome, EVENT).
    let mut fives: HashMap<&[u8], Vec<&Inv>> = HashMap::new();
    let mut by_event: HashMap<(&[u8], &Value), Vec<&Inv>> = HashMap::new();
    for b in invs.iter().filter(|b| b.side == Side::Five) {
        fives.entry(chroms[b.n]).or_default().push(b);
        if let Some(e) = &b.event {
            by_event.entry((chroms[b.n], e)).or_default().push(b);
        }
    }
    for v in fives.values_mut() {
        v.sort_by_key(|b| b.pos);
    }
    for a in invs.iter().filter(|a| a.side == Side::Three) {
        let chrom = chroms[a.n];
        let mark = |b: &Inv, out: &mut Vec<bool>| {
            out[a.n] = true;
            out[b.n] = true;
        };
        if let Some(v) = fives.get(chrom)
            && pair_bp >= 0
        {
            let from = v.partition_point(|b| b.pos < a.pos - pair_bp);
            for b in v[from..].iter().take_while(|b| b.pos <= a.pos + pair_bp) {
                if (a.end - b.end).abs() <= pair_bp {
                    mark(b, &mut out);
                }
            }
        }
        if let Some(e) = &a.event
            && let Some(v) = by_event.get(&(chrom, e))
        {
            for b in v {
                mark(b, &mut out);
            }
        }
    }
    out
}

/// What the filter did, for the summary line.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub paired: usize,
    pub single: usize,
    pub short: usize,
    pub long: usize,
}

/// Filters `converted` with Manta's own `raw` records into `out`.
pub fn run(raw: &str, converted: &str, b: Bounds, mut out: impl Write) -> Result<Counts, InvError> {
    let rec_err =
        |path: &str, line: u64, why: &str| InvError::Record(path.to_owned(), line, why.to_owned());
    let raw_vcf = read_vcf(raw)?;
    let mut raw_by_pair: HashMap<Vec<u8>, Vec<&Record>> = HashMap::new();
    for r in &raw_vcf.records {
        let id =
            r.1.get(2)
                .ok_or_else(|| rec_err(raw, r.0, "a record without an ID column"))?;
        if let Some(k) = pair_key(id) {
            raw_by_pair.entry(k.to_vec()).or_default().push(r);
        }
    }
    let Vcf {
        header,
        mut records,
    } = read_vcf(converted)?;

    let mut invs = Vec::new();
    for (n, (line, f)) in records.iter().enumerate() {
        if f.len() < 8 {
            return Err(rec_err(
                converted,
                *line,
                &format!("a record with {} columns, not 8 or more", f.len()),
            ));
        }
        let d = info_dict(&f[7]);
        if d.get(&b"SVTYPE"[..]) != Some(&Value::Str(b"INV".to_vec()))
            || !f[2].starts_with(b"MantaINV")
        {
            continue;
        }
        let side = if d.contains_key(&b"INV3"[..]) {
            Side::Three
        } else if d.contains_key(&b"INV5"[..]) {
            Side::Five
        } else {
            continue;
        };
        // int(True) is 1 in Python: a bare END flag reads as 1.
        let end = match d.get(&b"END"[..]) {
            Some(Value::Flag) => 1,
            Some(Value::Str(s)) => match py2_int(s) {
                Some(e) => e,
                None => continue,
            },
            None => continue,
        };
        let pos = py2_int(&f[1]).ok_or_else(|| rec_err(converted, *line, &bad_pos(&f[1])))?;
        invs.push(Inv {
            n,
            pos,
            end,
            side,
            event: d.get(&b"EVENT"[..]).cloned(),
        });
    }
    let chroms: Vec<&[u8]> = records.iter().map(|(_, f)| f[0].as_slice()).collect();
    let paired = paired(&invs, &chroms, b.pair_bp);

    let mut c = Counts::default();
    // (line, fields) of what comes out: the replacement breakends first,
    // then the converted records not replaced, in file order.
    let mut out_recs: Vec<(u64, Vec<Vec<u8>>)> = Vec::new();
    let mut replaced = vec![false; records.len()];
    for inv in &invs {
        let f = &mut records[inv.n].1;
        if paired[inv.n] {
            add_info(f, b"INV_JUNCTIONS=2");
            c.paired += 1;
            continue;
        }
        let span = inv.end - inv.pos;
        let mates = raw_by_pair
            .get(&f[2][b"MantaINV".len()..])
            .map_or(&[][..], Vec::as_slice);
        if span >= b.max_single_bp && mates.len() == 2 {
            for (line, m) in mates {
                if m.len() < 8 {
                    return Err(rec_err(
                        raw,
                        *line,
                        &format!("a record with {} columns, not 8 or more", m.len()),
                    ));
                }
                let mut m = m.clone();
                add_info(&mut m, b"INV_JUNCTIONS=1");
                out_recs.push((*line, m));
            }
            replaced[inv.n] = true;
            c.long += 1;
            continue;
        }
        add_info(f, b"INV_JUNCTIONS=1");
        c.single += 1;
        if span < b.short_bp {
            add_info(f, b"INVDUP_CANDIDATE");
            c.short += 1;
        }
    }
    let from_raw = out_recs.len();
    out_recs.extend(
        records
            .into_iter()
            .zip(replaced)
            .filter(|(_, r)| !r)
            .map(|(rec, _)| rec),
    );

    // The ##contig order: a repeated contig takes its last place, as a dict
    // built from the list does; a contig not listed goes after, in the
    // order first seen, ranked by the dict's size at that moment.
    let mut order: HashMap<Vec<u8>, usize> = HashMap::new();
    let contigs = header
        .iter()
        .filter_map(|h| h.strip_prefix(b"##contig=<ID="));
    for (i, rest) in contigs.enumerate() {
        let id = rest.split(|&c| c == b',').next().unwrap_or_default();
        let id = id
            .iter()
            .rposition(|&c| c != b'>')
            .map_or(&id[..0], |k| &id[..=k]);
        order.insert(id.to_vec(), i);
    }
    let mut keyed = Vec::with_capacity(out_recs.len());
    for (i, (line, f)) in out_recs.into_iter().enumerate() {
        if !order.contains_key(&f[0]) {
            let at = order.len();
            order.insert(f[0].clone(), at);
        }
        let path = if i < from_raw { raw } else { converted };
        let pos = py2_int(&f[1]).ok_or_else(|| rec_err(path, line, &bad_pos(&f[1])))?;
        keyed.push((order[&f[0]], pos, f));
    }
    keyed.sort_by_key(|k| (k.0, k.1));

    let new_info = [
        format!(
            "##INFO=<ID=INV_JUNCTIONS,Number=1,Type=Integer,Description=\"Junctions Manta found for this inversion: 2 when the reciprocal INV3 and INV5 junctions are both present (POS and END within {} bp, or one Manta EVENT), 1 when only one is; a single junction of {} bp or more stays a breakend pair\">",
            b.pair_bp, b.max_single_bp
        ),
        format!(
            "##INFO=<ID=INVDUP_CANDIDATE,Number=0,Type=Flag,Description=\"Single-junction inversion shorter than {} bp: more likely an inverted insertion or duplication than an inversion\">",
            b.short_bp
        ),
    ];
    let w = |e| InvError::Io("<stdout>".to_owned(), e);
    for h in &header {
        if h.starts_with(b"#CHROM") {
            for l in &new_info {
                writeln!(out, "{l}").map_err(w)?;
            }
        }
        out.write_all(h).map_err(w)?;
        out.write_all(b"\n").map_err(w)?;
    }
    for (_, _, f) in keyed {
        out.write_all(&f.join(&b'\t')).map_err(w)?;
        out.write_all(b"\n").map_err(w)?;
    }
    out.flush().map_err(w)?;
    Ok(c)
}

/// `og sv inv-filter`: exit 0, or 2 for a usage error or an input it
/// cannot read.
pub fn cli(argv: &[String]) -> u8 {
    let parse = || -> Result<(String, String, Bounds), ArgError> {
        let a = args::parse(&SPEC, argv)?;
        let d = Bounds::default();
        let b = Bounds {
            pair_bp: a.typed("--pair-bp", d.pair_bp, "int", args::int)?,
            max_single_bp: a.typed("--max-single-bp", d.max_single_bp, "int", args::int)?,
            short_bp: a.typed("--short-bp", d.short_bp, "int", args::int)?,
        };
        Ok((a.pos[0].clone(), a.pos[1].clone(), b))
    };
    let (raw, converted, b) = match parse() {
        Ok(p) => p,
        Err(e) => return args::report(&SPEC, &e),
    };
    let stdout = io::BufWriter::with_capacity(1 << 16, io::stdout().lock());
    match run(&raw, &converted, b, stdout) {
        Ok(c) => {
            eprintln!(
                "Inversion check: {} INV records with both junctions, {} single-junction INV records kept ({} shorter than {} bp marked INVDUP_CANDIDATE), {} single-junction inversions of {} bp or more left as their breakend pairs",
                c.paired, c.single, c.short, b.short_bp, c.long, b.max_single_bp
            );
            0
        }
        Err(e) => {
            eprintln!("ERROR: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &std::path::Path, name: &str, text: &str) -> String {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn inv(chrom: &str, pos: u64, end: u64, id: &str, side: u8, event: Option<&str>) -> String {
        let ev = event.map(|e| format!(";EVENT={e}")).unwrap_or_default();
        format!(
            "{chrom}\t{pos}\tMantaINV:{id}\tN\t<INV>\t.\tPASS\tEND={end};SVTYPE=INV;INV{side}{ev}\tGT\t0/1\n"
        )
    }

    fn bnd(chrom: &str, pos: u64, id: &str, mate: u32) -> String {
        format!("{chrom}\t{pos}\tMantaBND:{id}:{mate}\tN\tN]x]\t.\tPASS\tSVTYPE=BND\tGT\t0/1\n")
    }

    const HEAD: &str = "##fileformat=VCFv4.1\n##contig=<ID=chr1,length=1>\n##contig=<ID=chr2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

    fn filter(raw: &str, conv: &str, b: Bounds) -> (String, Counts) {
        let dir = std::env::temp_dir().join(format!("og-inv-{}-{}", std::process::id(), rand_id()));
        std::fs::create_dir_all(&dir).unwrap();
        let r = write(&dir, "raw.vcf", raw);
        let c = write(&dir, "conv.vcf", conv);
        let mut out = Vec::new();
        let counts = run(&r, &c, b, &mut out).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        (String::from_utf8(out).unwrap(), counts)
    }

    fn rand_id() -> u64 {
        use std::sync::atomic::{AtomicU64, Ordering};
        static N: AtomicU64 = AtomicU64::new(0);
        N.fetch_add(1, Ordering::Relaxed)
    }

    #[test]
    fn pairs_singles_long_and_short() {
        let raw = format!(
            "{HEAD}{}{}",
            bnd("chr1", 5_000_000, ":9:0:1:0:0:0", 0),
            bnd("chr1", 60_000_000, ":9:0:1:0:0:0", 1)
        );
        let conv = format!(
            "{HEAD}{}{}{}{}{}",
            inv("chr1", 100, 2000, ":1:0:1:0:0:0", 3, None),
            inv("chr1", 200, 2100, ":2:0:1:0:0:0", 5, None),
            inv("chr1", 5_000_000, 60_000_000, ":9:0:1:0:0:0", 3, None),
            inv("chr2", 10, 410, ":3:0:1:0:0:0", 5, None),
            inv("chr2", 900, 200_000, ":4:0:1:0:0:0", 3, None),
        );
        let (out, c) = filter(&raw, &conv, Bounds::default());
        assert_eq!(
            c,
            Counts {
                paired: 2,
                single: 2,
                short: 1,
                long: 1
            }
        );
        let recs: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(recs.len(), 6);
        assert!(recs[0].starts_with("chr1\t100\t") && recs[0].contains("INV3;INV_JUNCTIONS=2"));
        assert!(
            recs[1].starts_with("chr1\t200\t") && recs[1].ends_with("INV_JUNCTIONS=2\tGT\t0/1")
        );
        assert!(
            recs[2].starts_with("chr1\t5000000\tMantaBND")
                && recs[2].contains("SVTYPE=BND;INV_JUNCTIONS=1")
        );
        assert!(recs[3].starts_with("chr1\t60000000\tMantaBND"));
        assert!(recs[4].contains("INV5;INV_JUNCTIONS=1;INVDUP_CANDIDATE"));
        assert!(recs[5].contains("INV3;INV_JUNCTIONS=1\t"));
        assert!(out.contains("##INFO=<ID=INVDUP_CANDIDATE,") && out.contains("#CHROM"));
    }

    #[test]
    fn red_first_the_pair_window_is_inclusive() {
        // 1000 bp apart in POS and END: paired at --pair-bp 1000, not 999.
        let conv = format!(
            "{HEAD}{}{}",
            inv("chr1", 100, 5000, ":1:0:1:0:0:0", 3, None),
            inv("chr1", 1100, 6000, ":2:0:1:0:0:0", 5, None),
        );
        let (_, c) = filter(HEAD, &conv, Bounds::default());
        assert_eq!(c.paired, 2);
        let b = Bounds {
            pair_bp: 999,
            ..Bounds::default()
        };
        let (_, c) = filter(HEAD, &conv, b);
        assert_eq!(c.paired, 0);
    }

    #[test]
    fn one_event_pairs_far_junctions() {
        let conv = format!(
            "{HEAD}{}{}",
            inv("chr1", 100, 5000, ":1:0:1:0:0:0", 3, Some("E1")),
            inv("chr1", 90_000, 95_000, ":2:0:1:0:0:0", 5, Some("E1")),
        );
        let (_, c) = filter(HEAD, &conv, Bounds::default());
        assert_eq!(c.paired, 2);
    }

    #[test]
    fn a_contig_not_in_the_header_goes_last() {
        let conv = format!(
            "{HEAD}chrZ\t5\t.\tN\t<DEL>\t.\tPASS\t.\nchr2\t9\t.\tN\t<DEL>\t.\tPASS\t.\nchr1\t7\t.\tN\t<DEL>\t.\tPASS\t.\n"
        );
        let (out, _) = filter(HEAD, &conv, Bounds::default());
        let chroms: Vec<&str> = out
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| &l[..4])
            .collect();
        assert_eq!(chroms, ["chr1", "chr2", "chrZ"]);
    }
}
