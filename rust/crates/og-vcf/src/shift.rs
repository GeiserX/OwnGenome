//! `og setup chrm-shift` and `og vcf shift-back`: the shifted-reference
//! pass for the chrM control region.
//!
//! chrM is a circle, but GRCh38 writes it as a line from m.1 to m.16569 (the
//! rCRS). A read across the origin, where m.16569 joins m.1, aligns to
//! neither end in full, so calls in the control region around the origin
//! come from fewer reads. Broad's mitochondria pipeline calls the control
//! region a second time on chrM rotated by 8000 bp, where the origin lies in
//! the middle of the sequence. Two tools do that here:
//!
//! - `og setup chrm-shift` writes the rotated chrM (its base 1 is m.8001),
//!   its `.fai`, and the control region in its positions for Mutect2's `-L`.
//! - `og vcf shift-back` maps the calls of the rotated pass back to chrM's
//!   own positions, `(pos - 1 + 8000) mod 16569 + 1`, so no chain file is
//!   needed, and merges them with the calls of the linear pass: the control
//!   region (m.16025-16569 and m.1-575, Broad's split) from the rotated pass,
//!   the rest from the linear one. The merged VCF says so in its
//!   `##og_shifted_pass=` header line.
//!
//! On the rotated chrM the control region is one stretch, m.16025 at 8025 to
//! m.575 at 9144, with the origin between 8569 (m.16569) and 8570 (m.1). A
//! call of the rotated pass whose reference bases cross that point cannot be
//! written on the linear chrM: it is dropped and counted on stderr.

use og_io::AtomicFile;
use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

/// The contig of the pass.
pub const CONTIG: &str = "chrM";
/// The length of GRCh38's chrM, the rCRS.
pub const RCRS_LEN: u64 = 16569;
/// How far the rotated chrM is turned: its base 1 is m.8001.
pub const SHIFT: u64 = 8000;
/// The control region of the rCRS as Broad's pipeline splits chrM: from
/// m.16025 across the origin to m.575.
pub const CONTROL: (u64, u64) = (16025, 575);
/// The header key of the line a merged VCF carries.
pub const MARK: &str = "##og_shifted_pass=";
/// Bases per line of the rotated FASTA.
const LINE: usize = 60;

/// The position on the rotated sequence of position `pos` (1-based, at most
/// `len`) of the original.
pub fn to_rotated(pos: u64, len: u64, shift: u64) -> u64 {
    (pos - 1 + len - shift % len) % len + 1
}

/// The position on the original of position `pos` (1-based, at most `len`)
/// of the rotated sequence.
pub fn from_rotated(pos: u64, len: u64, shift: u64) -> u64 {
    (pos - 1 + shift % len) % len + 1
}

/// `seq` rotated by `shift`: base `i` of the result is base
/// `(i + shift) mod len` of `seq`, so position `p` of the result is position
/// [`from_rotated`]`(p)` of `seq`.
pub fn rotate(seq: &[u8], shift: u64) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }
    let k = usize::try_from(shift % seq.len() as u64).expect("below the length");
    let mut out = Vec::with_capacity(seq.len());
    out.extend_from_slice(&seq[k..]);
    out.extend_from_slice(&seq[..k]);
    out
}

/// Whether rCRS position `pos` is in the control region.
pub fn in_control(pos: u64) -> bool {
    pos >= CONTROL.0 || pos <= CONTROL.1
}

/// The control region on the rotated chrM, one stretch: 8025-9144.
pub fn control_rotated() -> (u64, u64) {
    (
        to_rotated(CONTROL.0, RCRS_LEN, SHIFT),
        to_rotated(CONTROL.1, RCRS_LEN, SHIFT),
    )
}

/// The value after `--key`, as `--key VALUE` or `--key=VALUE`.
fn take(
    a: &str,
    it: &mut std::slice::Iter<'_, String>,
    known: &[&str],
) -> Result<(String, String), String> {
    let (key, inline) = match a.split_once('=') {
        Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
        _ => (a, None),
    };
    if !known.contains(&key) {
        return Err(format!("unrecognized arguments: {a}"));
    }
    let v = inline
        .or_else(|| it.next().cloned())
        .ok_or_else(|| format!("argument {key}: expected one argument"))?;
    Ok((key.to_owned(), v))
}

// ---------------------------------------------------------------------------
// og setup chrm-shift

const SETUP_USAGE: &str =
    "usage: og setup chrm-shift --fasta REF.fa --out OUT.fasta [--intervals OUT.intervals]\n";

/// One contig's line of a `.fai`.
#[derive(Debug, PartialEq)]
struct Fai {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

fn fai_entry(fai: &Path, name: &str) -> Result<Fai, String> {
    let text = std::fs::read_to_string(fai).map_err(|e| format!("{}: {e}", fai.display()))?;
    for line in text.lines() {
        let c: Vec<&str> = line.split('\t').collect();
        if c.first() != Some(&name) {
            continue;
        }
        let num = |i: usize| {
            c.get(i).and_then(|v| v.parse::<u64>().ok()).ok_or_else(|| {
                format!("{}: the {name} line has no column {}", fai.display(), i + 1)
            })
        };
        let e = Fai {
            len: num(1)?,
            offset: num(2)?,
            line_bases: num(3)?,
            line_width: num(4)?,
        };
        if e.line_bases == 0 || e.line_width < e.line_bases {
            return Err(format!(
                "{}: the {name} line is not a FASTA index line",
                fai.display()
            ));
        }
        return Ok(e);
    }
    Err(format!("{}: no {name} in the index", fai.display()))
}

/// The bases of the contig `e` points at, read from its offset.
fn read_contig(fasta: &Path, e: &Fai) -> Result<Vec<u8>, String> {
    let fail = |m: String| format!("{}: {m}", fasta.display());
    let mut f = File::open(fasta).map_err(|x| fail(x.to_string()))?;
    f.seek(SeekFrom::Start(e.offset))
        .map_err(|x| fail(x.to_string()))?;
    let bytes = e.len / e.line_bases * e.line_width + e.len % e.line_bases;
    let mut raw = Vec::with_capacity(usize::try_from(bytes).unwrap_or(0));
    f.take(bytes)
        .read_to_end(&mut raw)
        .map_err(|x| fail(x.to_string()))?;
    let seq: Vec<u8> = raw
        .into_iter()
        .filter(|b| !matches!(b, b'\n' | b'\r'))
        .collect();
    if seq.len() as u64 != e.len || seq.iter().any(|b| !b.is_ascii_alphabetic()) {
        return Err(fail(format!(
            "the {CONTIG} sequence at offset {} is not the {} bases its .fai names: rebuild the .fai",
            e.offset, e.len
        )));
    }
    Ok(seq)
}

fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut f = AtomicFile::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    f.write_all(bytes)
        .and_then(|()| f.commit())
        .map_err(|e| format!("{}: {e}", path.display()))
}

/// `path` with `suffix` added to its name.
fn plus(path: &Path, suffix: &str) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(suffix);
    PathBuf::from(p)
}

/// Writes the rotated chrM of `fasta` to `out` with its `.fai`, and the
/// control region in its positions to `intervals`.
fn chrm_shift(fasta: &Path, out: &Path, intervals: Option<&Path>) -> Result<(), String> {
    let e = fai_entry(&plus(fasta, ".fai"), CONTIG)?;
    if e.len != RCRS_LEN {
        return Err(format!(
            "{}: {CONTIG} is {} bp, not the {RCRS_LEN} bp of the rCRS: the control region is \
             defined on the rCRS, so this reference has no shifted pass",
            fasta.display(),
            e.len
        ));
    }
    let seq = rotate(&read_contig(fasta, &e)?, SHIFT);
    let head = format!(">{CONTIG}\n");
    let mut text = head.clone().into_bytes();
    for line in seq.chunks(LINE) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
    let fai = format!(
        "{CONTIG}\t{}\t{}\t{LINE}\t{}\n",
        seq.len(),
        head.len(),
        LINE + 1
    );
    write_atomic(out, &text)?;
    write_atomic(&plus(out, ".fai"), fai.as_bytes())?;
    if let Some(p) = intervals {
        let (a, b) = control_rotated();
        write_atomic(p, format!("{CONTIG}:{a}-{b}\n").as_bytes())?;
    }
    Ok(())
}

/// `og setup chrm-shift --fasta REF.fa --out OUT.fasta [--intervals FILE]`:
/// writes chrM of `REF.fa` (read through `REF.fa.fai`) rotated by 8000 bp,
/// its `.fai`, and the control region on it (`chrM:8025-9144`) for
/// Mutect2's `-L`. Exit 0 done, 1 an error, 2 usage.
pub fn chrm_shift_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut usage = |m: &str| {
        let _ = writeln!(err, "{SETUP_USAGE}og setup chrm-shift: error: {m}");
        2
    };
    let (mut fasta, mut dest, mut intervals) = (None, None, None);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if a == "-h" || a == "--help" {
            let _ = write!(out, "{SETUP_USAGE}");
            return 0;
        }
        match take(a, &mut it, &["--fasta", "--out", "--intervals"]) {
            Ok((k, v)) => match k.as_str() {
                "--fasta" => fasta = Some(PathBuf::from(v)),
                "--out" => dest = Some(PathBuf::from(v)),
                _ => intervals = Some(PathBuf::from(v)),
            },
            Err(m) => return usage(&m),
        }
    }
    let (Some(fasta), Some(dest)) = (fasta, dest) else {
        return usage("the following arguments are required: --fasta, --out");
    };
    match chrm_shift(&fasta, &dest, intervals.as_deref()) {
        Ok(()) => 0,
        Err(m) => {
            let _ = writeln!(err, "og setup chrm-shift: {m}");
            1
        }
    }
}

// ---------------------------------------------------------------------------
// og vcf shift-back

const BACK_USAGE: &str =
    "usage: og vcf shift-back --linear LINEAR.vcf.gz --shifted SHIFTED.vcf.gz --out OUT.vcf.gz\n";

/// A VCF read whole: the `##` lines, the `#CHROM` line, the records.
struct Vcf {
    meta: Vec<Vec<u8>>,
    columns: Vec<u8>,
    records: Vec<Vec<u8>>,
}

fn read_vcf(path: &Path) -> Result<Vcf, String> {
    let fail = |e: io::Error| format!("{}: {e}", path.display());
    let mut r = og_io::open(path).map_err(fail)?;
    let mut v = Vcf {
        meta: Vec::new(),
        columns: Vec::new(),
        records: Vec::new(),
    };
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if r.read_until(b'\n', &mut buf).map_err(fail)? == 0 {
            break;
        }
        let mut line = buf.strip_suffix(b"\n").unwrap_or(&buf);
        line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.starts_with(b"##") {
            v.meta.push(line.to_vec());
        } else if line.starts_with(b"#") {
            v.columns = line.to_vec();
        } else if !line.is_empty() {
            v.records.push(line.to_vec());
        }
    }
    if v.columns.is_empty() {
        return Err(format!("{}: no #CHROM line: not a VCF", path.display()));
    }
    Ok(v)
}

/// The length a `##contig` line gives `CONTIG`, if there is one.
fn contig_len(meta: &[Vec<u8>]) -> Option<u64> {
    let want = format!("##contig=<ID={CONTIG},");
    meta.iter().find_map(|l| {
        let rest = l.strip_prefix(want.as_bytes())?;
        let rest = std::str::from_utf8(rest).ok()?.trim_end_matches('>');
        rest.split(',')
            .find_map(|kv| kv.strip_prefix("length="))
            .and_then(|n| n.parse().ok())
    })
}

/// The ID of a `##FILTER`, `##INFO` or `##FORMAT` line, with its kind.
fn def_id(line: &[u8]) -> Option<(&[u8], &[u8])> {
    let s = line.strip_prefix(b"##")?;
    let eq = s.iter().position(|b| *b == b'=')?;
    let (kind, rest) = (&s[..eq], &s[eq + 1..]);
    if !matches!(kind, b"FILTER" | b"INFO" | b"FORMAT") {
        return None;
    }
    let rest = rest.strip_prefix(b"<ID=")?;
    let end = rest.iter().position(|b| matches!(b, b',' | b'>'))?;
    Some((kind, &rest[..end]))
}

/// The position of a record of `path`, checked to be a chrM position.
fn record_pos(path: &Path, cols: &[&[u8]]) -> Result<u64, String> {
    if cols.len() < 8 {
        return Err(format!(
            "{}: a record has {} columns, not 8 or more",
            path.display(),
            cols.len()
        ));
    }
    if cols[0] != CONTIG.as_bytes() {
        return Err(format!(
            "{}: a record on {}: only {CONTIG} calls are merged",
            path.display(),
            String::from_utf8_lossy(cols[0])
        ));
    }
    std::str::from_utf8(cols[1])
        .ok()
        .and_then(|p| p.parse::<u64>().ok())
        .filter(|p| (1..=RCRS_LEN).contains(p))
        .ok_or_else(|| {
            format!(
                "{}: POS {} is not a position of the {RCRS_LEN} bp {CONTIG}",
                path.display(),
                String::from_utf8_lossy(cols[1])
            )
        })
}

/// A sample value of FORMAT key `key` mapped back from the rotated chrM:
/// `PS` is a position, `PID` starts with one (`POS_REF_ALT`).
fn map_format(key: &[u8], value: &[u8]) -> Vec<u8> {
    let back = |p: &[u8]| -> Option<u64> {
        std::str::from_utf8(p)
            .ok()?
            .parse::<u64>()
            .ok()
            .filter(|p| (1..=RCRS_LEN).contains(p))
            .map(|p| from_rotated(p, RCRS_LEN, SHIFT))
    };
    match key {
        b"PS" => back(value).map_or_else(|| value.to_vec(), |p| p.to_string().into_bytes()),
        b"PID" => {
            let cut = value.iter().position(|b| *b == b'_').unwrap_or(value.len());
            match back(&value[..cut]) {
                Some(p) => [p.to_string().as_bytes(), &value[cut..]].concat(),
                None => value.to_vec(),
            }
        }
        _ => value.to_vec(),
    }
}

/// What happened to a record of the rotated pass.
enum Back {
    /// In the control region: its position on chrM and the rewritten line.
    Kept(u64, Vec<u8>),
    /// Outside the control region: the linear pass calls it.
    Outside,
    /// Its reference bases cross the origin.
    Crosses,
}

fn shift_record(path: &Path, line: &[u8]) -> Result<Back, String> {
    let cols: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
    let s = record_pos(path, &cols)?;
    let last = s + cols[3].len().max(1) as u64 - 1;
    let origin = RCRS_LEN - SHIFT; // m.16569 on the rotated chrM
    if s <= origin && last > origin {
        return Ok(Back::Crosses);
    }
    let pos = from_rotated(s, RCRS_LEN, SHIFT);
    if !in_control(pos) {
        return Ok(Back::Outside);
    }
    let keys: Vec<&[u8]> = cols
        .get(8)
        .map_or_else(Vec::new, |f| f.split(|b| *b == b':').collect());
    let mut out = Vec::with_capacity(line.len() + 2);
    for (i, c) in cols.iter().enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        match i {
            1 => out.extend_from_slice(pos.to_string().as_bytes()),
            _ if i > 8 && (keys.contains(&&b"PS"[..]) || keys.contains(&&b"PID"[..])) => {
                for (j, v) in c.split(|b| *b == b':').enumerate() {
                    if j > 0 {
                        out.push(b':');
                    }
                    out.extend_from_slice(&map_format(keys.get(j).copied().unwrap_or(b""), v));
                }
            }
            _ => out.extend_from_slice(c),
        }
    }
    Ok(Back::Kept(pos, out))
}

/// How many records of each kind went into a merge.
#[derive(Debug, Default, PartialEq)]
struct Counts {
    shifted: usize,
    linear: usize,
    crossing: usize,
}

/// The merged VCF of the linear and the rotated pass, as text, and its counts.
fn merge(linear_path: &Path, shifted_path: &Path) -> Result<(Vec<u8>, Counts), String> {
    let linear = read_vcf(linear_path)?;
    let shifted = read_vcf(shifted_path)?;
    for (p, v) in [(linear_path, &linear), (shifted_path, &shifted)] {
        if let Some(n) = contig_len(&v.meta).filter(|n| *n != RCRS_LEN) {
            return Err(format!(
                "{}: {CONTIG} is {n} bp, not the {RCRS_LEN} bp of the rCRS",
                p.display()
            ));
        }
    }
    if linear.columns != shifted.columns {
        return Err(format!(
            "{} and {} name different samples: both passes must call the same reads",
            linear_path.display(),
            shifted_path.display()
        ));
    }
    let mut counts = Counts::default();
    let mut records: Vec<(u64, Vec<u8>)> = Vec::new();
    for line in &linear.records {
        let cols: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
        let pos = record_pos(linear_path, &cols)?;
        if !in_control(pos) {
            records.push((pos, line.clone()));
            counts.linear += 1;
        }
    }
    for line in &shifted.records {
        match shift_record(shifted_path, line)? {
            Back::Kept(pos, l) => {
                records.push((pos, l));
                counts.shifted += 1;
            }
            Back::Outside => {}
            Back::Crosses => counts.crossing += 1,
        }
    }
    records.sort_by_key(|(p, _)| *p);
    let known: Vec<(&[u8], &[u8])> = linear.meta.iter().filter_map(|l| def_id(l)).collect();
    let mut text = Vec::new();
    for l in &linear.meta {
        text.extend_from_slice(l);
        text.push(b'\n');
    }
    for l in &shifted.meta {
        if def_id(l).is_some_and(|d| !known.contains(&d)) {
            text.extend_from_slice(l);
            text.push(b'\n');
        }
    }
    let (a, b) = CONTROL;
    text.extend_from_slice(
        format!(
            "{MARK}{CONTIG}:{a}-{RCRS_LEN} and {CONTIG}:1-{b} from Mutect2 on {CONTIG} rotated by \
             {SHIFT} bp, mapped back by og vcf shift-back; the rest from the linear pass\n"
        )
        .as_bytes(),
    );
    text.extend_from_slice(&linear.columns);
    text.push(b'\n');
    for (_, l) in &records {
        text.extend_from_slice(l);
        text.push(b'\n');
    }
    Ok((text, counts))
}

/// Writes `text` to `out` as BGZF, and its `.tbi` beside it.
fn write_bgzf_indexed(out: &Path, text: &[u8]) -> Result<(), String> {
    let fail = |e: io::Error| format!("{}: {e}", out.display());
    let mut w = og_io::bgzf::Writer::new(Vec::new(), NonZero::<usize>::MIN).map_err(fail)?;
    w.write_all(text).map_err(fail)?;
    let bytes = w.finish().map_err(fail)?;
    let index = og_io::tabix_write::build(io::Cursor::new(&bytes)).map_err(fail)?;
    write_atomic(out, &bytes)?;
    let tbi = plus(out, ".tbi");
    index
        .write_tbi(&tbi, true)
        .map_err(|e| format!("{}: {e}", tbi.display()))
}

/// `og vcf shift-back --linear L.vcf.gz --shifted S.vcf.gz --out OUT.vcf.gz`:
/// the control-region calls of the rotated pass `S`, mapped back to chrM,
/// merged with the calls of the linear pass `L` outside the control region,
/// into `OUT` (BGZF) and `OUT.tbi`. Exit 0 done, 1 an error, 2 usage.
pub fn shift_back_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut usage = |m: &str| {
        let _ = writeln!(err, "{BACK_USAGE}og vcf shift-back: error: {m}");
        2
    };
    let (mut linear, mut shifted, mut dest) = (None, None, None);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if a == "-h" || a == "--help" {
            let _ = write!(out, "{BACK_USAGE}");
            return 0;
        }
        match take(a, &mut it, &["--linear", "--shifted", "--out"]) {
            Ok((k, v)) => match k.as_str() {
                "--linear" => linear = Some(PathBuf::from(v)),
                "--shifted" => shifted = Some(PathBuf::from(v)),
                _ => dest = Some(PathBuf::from(v)),
            },
            Err(m) => return usage(&m),
        }
    }
    let (Some(linear), Some(shifted), Some(dest)) = (linear, shifted, dest) else {
        return usage("the following arguments are required: --linear, --shifted, --out");
    };
    let done = merge(&linear, &shifted).and_then(|(text, c)| {
        write_bgzf_indexed(&dest, &text)?;
        Ok(c)
    });
    match done {
        Ok(c) => {
            let _ = writeln!(
                err,
                "og vcf shift-back: {} control-region calls from the rotated pass, {} calls of the \
                 linear pass outside it, {} rotated calls across the origin dropped",
                c.shifted, c.linear, c.crossing
            );
            0
        }
        Err(m) => {
            let _ = writeln!(err, "og vcf shift-back: {m}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A fixed, uneven chrM-like sequence of `n` bases.
    fn bases(n: usize) -> Vec<u8> {
        let mut x: u32 = 0x9e37_79b9;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 17;
                x ^= x << 5;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn every_chrm_position_goes_round_and_back() {
        let mut seen = vec![false; RCRS_LEN as usize + 1];
        for p in 1..=RCRS_LEN {
            let r = to_rotated(p, RCRS_LEN, SHIFT);
            assert!((1..=RCRS_LEN).contains(&r), "m.{p} -> {r}");
            assert_eq!(from_rotated(r, RCRS_LEN, SHIFT), p, "m.{p} -> {r} -> back");
            assert_eq!(
                to_rotated(from_rotated(p, RCRS_LEN, SHIFT), RCRS_LEN, SHIFT),
                p
            );
            assert!(!seen[r as usize], "two positions rotate to {r}");
            seen[r as usize] = true;
        }
        // The plan's formula and the ends of the origin.
        for s in 1..=RCRS_LEN {
            assert_eq!(from_rotated(s, RCRS_LEN, SHIFT), (s - 1 + 8000) % 16569 + 1);
        }
        assert_eq!(from_rotated(1, RCRS_LEN, SHIFT), 8001);
        assert_eq!(from_rotated(8569, RCRS_LEN, SHIFT), 16569);
        assert_eq!(from_rotated(8570, RCRS_LEN, SHIFT), 1);
    }

    #[test]
    fn every_position_goes_round_and_back_for_any_length_and_shift() {
        for len in 1..=40u64 {
            for shift in 0..=2 * len {
                for p in 1..=len {
                    let r = to_rotated(p, len, shift);
                    assert!((1..=len).contains(&r));
                    assert_eq!(
                        from_rotated(r, len, shift),
                        p,
                        "len {len} shift {shift} pos {p}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_rotated_sequence_holds_each_base_at_its_rotated_position() {
        let seq = bases(RCRS_LEN as usize);
        let rot = rotate(&seq, SHIFT);
        assert_eq!(rot.len(), seq.len());
        for s in 1..=RCRS_LEN {
            let p = from_rotated(s, RCRS_LEN, SHIFT);
            assert_eq!(
                rot[s as usize - 1],
                seq[p as usize - 1],
                "rotated {s} is m.{p}"
            );
        }
        assert_eq!(
            rotate(&rot, RCRS_LEN - SHIFT),
            seq,
            "turning the rest of the way gives chrM back"
        );
        assert_eq!(rotate(&seq, 0), seq);
        assert!(rotate(b"", SHIFT).is_empty());
    }

    #[test]
    fn the_control_region_is_one_stretch_on_the_rotated_chrm() {
        assert_eq!(control_rotated(), (8025, 9144));
        for p in 1..=RCRS_LEN {
            let r = to_rotated(p, RCRS_LEN, SHIFT);
            assert_eq!(in_control(p), (8025..=9144).contains(&r), "m.{p} at {r}");
        }
        assert!(in_control(16025) && in_control(16569) && in_control(1) && in_control(575));
        assert!(!in_control(16024) && !in_control(576));
    }

    /// A scratch directory removed on drop.
    struct Tmp(PathBuf);
    impl Tmp {
        fn new(name: &str) -> Tmp {
            let d = std::env::temp_dir().join(format!("og-shift-{name}-{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&d);
            std::fs::create_dir_all(&d).unwrap();
            Tmp(d)
        }
        fn put(&self, name: &str, text: &[u8]) -> PathBuf {
            let p = self.0.join(name);
            std::fs::write(&p, text).unwrap();
            p
        }
    }
    impl Drop for Tmp {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn run(f: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8, args: &[&str]) -> (u8, String) {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let rc = f(&args, &mut out, &mut err);
        (rc, String::from_utf8(err).unwrap())
    }

    /// A FASTA with a contig before chrM, lines of `width`, and its `.fai`.
    fn reference(t: &Tmp, chrm: &[u8], width: usize) -> PathBuf {
        let mut text = b">chr1 first\nNNNNACGT\nAC\n".to_vec();
        let mut fai = String::from("chr1\t10\t12\t8\t9\n");
        let head = b">chrM description\n";
        fai.push_str(&format!(
            "chrM\t{}\t{}\t{width}\t{}\n",
            chrm.len(),
            text.len() + head.len(),
            width + 1
        ));
        text.extend_from_slice(head);
        for l in chrm.chunks(width) {
            text.extend_from_slice(l);
            text.push(b'\n');
        }
        t.put("ref.fa.fai", fai.as_bytes());
        t.put("ref.fa", &text)
    }

    #[test]
    fn chrm_shift_writes_the_rotated_chrm_its_index_and_the_control_interval() {
        let t = Tmp::new("setup");
        let chrm = bases(RCRS_LEN as usize);
        let fa = reference(&t, &chrm, 70);
        let (out, iv) = (
            t.0.join("chrM_shifted.fasta"),
            t.0.join("control.intervals"),
        );
        let (rc, err) = run(
            chrm_shift_cli,
            &[
                "--fasta",
                fa.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
                "--intervals",
                iv.to_str().unwrap(),
            ],
        );
        assert_eq!((rc, err.as_str()), (0, ""));
        let text = std::fs::read(&out).unwrap();
        let mut want = b">chrM\n".to_vec();
        for l in rotate(&chrm, SHIFT).chunks(60) {
            want.extend_from_slice(l);
            want.push(b'\n');
        }
        assert_eq!(text, want);
        assert_eq!(
            std::fs::read_to_string(plus(&out, ".fai")).unwrap(),
            "chrM\t16569\t6\t60\t61\n"
        );
        assert_eq!(std::fs::read_to_string(&iv).unwrap(), "chrM:8025-9144\n");
        // The written .fai reads the rotated chrM back.
        let e = fai_entry(&plus(&out, ".fai"), CONTIG).unwrap();
        assert_eq!(read_contig(&out, &e).unwrap(), rotate(&chrm, SHIFT));
    }

    #[test]
    fn chrm_shift_refuses_a_chrm_that_is_not_the_rcrs_or_absent() {
        let t = Tmp::new("setup-bad");
        let fa = reference(&t, &bases(16571), 60);
        let out = t.0.join("o.fasta");
        let (rc, err) = run(
            chrm_shift_cli,
            &[
                "--fasta",
                fa.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 1);
        assert!(
            err.contains("chrM is 16571 bp, not the 16569 bp of the rCRS"),
            "{err}"
        );
        assert!(!out.exists());
        t.put("ref.fa.fai", b"chr1\t10\t12\t8\t9\n");
        let (rc, err) = run(
            chrm_shift_cli,
            &[
                "--fasta",
                fa.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 1);
        assert!(err.contains("no chrM in the index"), "{err}");
        // A stale index that points into the wrong bytes.
        t.put("ref.fa.fai", b"chrM\t16569\t3\t60\t61\n");
        let (rc, err) = run(
            chrm_shift_cli,
            &[
                "--fasta",
                fa.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 1);
        assert!(err.contains("rebuild the .fai"), "{err}");
        let (rc, err) = run(chrm_shift_cli, &["--fasta", fa.to_str().unwrap()]);
        assert_eq!(rc, 2);
        assert!(err.contains("required: --fasta, --out"), "{err}");
    }

    const META: &str = "##fileformat=VCFv4.2\n\
        ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
        ##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fraction\">\n\
        ##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Phase set id\">\n\
        ##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n\
        ##contig=<ID=chrM,length=16569>\n";
    const COLS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

    fn rec(pos: u64, r: &str, a: &str, filter: &str, sample: &str) -> String {
        format!("chrM\t{pos}\t.\t{r}\t{a}\t.\t{filter}\tDP=900\tGT:AF:PID:PS\t{sample}\n")
    }

    fn rot(p: u64) -> u64 {
        to_rotated(p, RCRS_LEN, SHIFT)
    }

    fn records(path: &Path) -> Vec<String> {
        let v = read_vcf(path).unwrap();
        v.records
            .into_iter()
            .map(|l| String::from_utf8(l).unwrap())
            .collect()
    }

    #[test]
    fn shift_back_takes_the_control_region_from_the_rotated_pass() {
        let t = Tmp::new("back");
        // The linear pass: m.73 and m.16519 in the control region (replaced
        // by the rotated pass), m.3243 and m.576 outside it (kept). It has no
        // call at m.3: the reads across the origin hid it.
        let linear = t.put(
            "linear.vcf",
            format!(
                "{META}{COLS}{}{}{}{}",
                rec(73, "A", "G", "PASS", "0/1:0.6:.:."),
                rec(576, "C", "T", "PASS", "0/1:0.5:.:."),
                rec(3243, "A", "G", "PASS", "0/1:0.3:.:."),
                rec(16519, "T", "C", "weak_evidence", "0/1:0.2:.:."),
            )
            .as_bytes(),
        );
        // The rotated pass, in its own positions and order: m.16519, m.16569
        // (its last base), a deletion across the origin, m.3 with a phase
        // set of m.16569, m.73, and m.3243 (outside, dropped).
        let shifted = t.put(
            "shifted.vcf",
            format!(
                "{META}{COLS}{}{}{}{}{}{}",
                rec(rot(3243), "A", "G", "PASS", "0/1:0.31:.:."),
                rec(rot(16519), "T", "C", "PASS", "1/1:0.99:.:."),
                rec(rot(16569), "C", "T", "PASS", "0|1:0.4:8569_C_T:8569"),
                rec(rot(16569), "CG", "C", "PASS", "0/1:0.1:.:."),
                rec(rot(3), "T", "C", "PASS", "0|1:0.4:8569_C_T:8569"),
                rec(rot(73), "A", "G", "PASS", "1/1:0.98:.:."),
            )
            .as_bytes(),
        );
        let out = t.0.join("merged.vcf.gz");
        let (rc, err) = run(
            shift_back_cli,
            &[
                "--linear",
                linear.to_str().unwrap(),
                "--shifted",
                shifted.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 0, "{err}");
        assert_eq!(
            err,
            "og vcf shift-back: 4 control-region calls from the rotated pass, 2 calls of the linear \
             pass outside it, 1 rotated calls across the origin dropped\n"
        );
        assert_eq!(
            records(&out),
            [
                rec(3, "T", "C", "PASS", "0|1:0.4:16569_C_T:16569"),
                rec(73, "A", "G", "PASS", "1/1:0.98:.:."),
                rec(576, "C", "T", "PASS", "0/1:0.5:.:."),
                rec(3243, "A", "G", "PASS", "0/1:0.3:.:."),
                rec(16519, "T", "C", "PASS", "1/1:0.99:.:."),
                rec(16569, "C", "T", "PASS", "0|1:0.4:16569_C_T:16569"),
            ]
            .map(|r| r.trim_end().to_owned())
        );
        let v = read_vcf(&out).unwrap();
        let meta: Vec<String> = v
            .meta
            .iter()
            .map(|l| String::from_utf8_lossy(l).into_owned())
            .collect();
        assert_eq!(
            meta[..6].join("\n") + "\n",
            META,
            "the linear pass's header comes first"
        );
        assert!(meta[6].starts_with(MARK), "{meta:?}");
        assert_eq!(meta.len(), 7);
        // BGZF with a tabix index that finds the records.
        assert_eq!(&std::fs::read(&out).unwrap()[..4], &[0x1f, 0x8b, 8, 4]);
        let mut r = og_io::tabix::Reader::open(&out).unwrap();
        let hits: Vec<Vec<u8>> = r.query("chrM", 0, 100).collect::<io::Result<_>>().unwrap();
        assert_eq!(hits.len(), 2, "m.3 and m.73");
    }

    #[test]
    fn shift_back_keeps_definitions_only_the_rotated_pass_has_and_refuses_mismatches() {
        let t = Tmp::new("back-head");
        let linear = t.put("l.vcf", format!("{META}{COLS}").as_bytes());
        let extra = "##FILTER=<ID=possible_numt,Description=\"NuMT\">";
        let shifted = t.put(
            "s.vcf",
            format!(
                "{META}{extra}\n{COLS}{}",
                rec(rot(100), "A", "G", "possible_numt", "0/1:0.1:.:.")
            )
            .as_bytes(),
        );
        let out = t.0.join("m.vcf.gz");
        let args = |l: &Path, s: &Path| {
            [
                "--linear".to_owned(),
                l.display().to_string(),
                "--shifted".to_owned(),
                s.display().to_string(),
                format!("--out={}", out.display()),
            ]
        };
        let (mut o, mut e) = (Vec::new(), Vec::new());
        assert_eq!(shift_back_cli(&args(&linear, &shifted), &mut o, &mut e), 0);
        let meta: Vec<Vec<u8>> = read_vcf(&out).unwrap().meta;
        assert!(
            meta.iter().any(|l| l == extra.as_bytes()),
            "the possible_numt definition is kept"
        );
        assert_eq!(
            records(&out),
            [rec(100, "A", "G", "possible_numt", "0/1:0.1:.:.")
                .trim_end()
                .to_owned()]
        );

        let other = t.put(
            "o.vcf",
            format!("{META}{}", COLS.replace("\tS\n", "\tT\n")).as_bytes(),
        );
        let (rc, err) = run(
            shift_back_cli,
            &[
                "--linear",
                linear.to_str().unwrap(),
                "--shifted",
                other.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 1);
        assert!(err.contains("name different samples"), "{err}");
        let hg19 = t.put(
            "h.vcf",
            format!("{}{COLS}", META.replace("16569", "16571")).as_bytes(),
        );
        let (rc, err) = run(
            shift_back_cli,
            &[
                "--linear",
                hg19.to_str().unwrap(),
                "--shifted",
                shifted.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 1);
        assert!(err.contains("chrM is 16571 bp"), "{err}");
        let chr1 = t.put(
            "c.vcf",
            format!(
                "{META}{COLS}{}",
                rec(5, "A", "G", "PASS", "0/1:0.1:.:.").replace("chrM", "chr1")
            )
            .as_bytes(),
        );
        let (rc, err) = run(
            shift_back_cli,
            &[
                "--linear",
                chr1.to_str().unwrap(),
                "--shifted",
                shifted.to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
            ],
        );
        assert_eq!(rc, 1);
        assert!(err.contains("a record on chr1"), "{err}");
        let (rc, err) = run(shift_back_cli, &["--linear", "x", "--bogus", "y"]);
        assert_eq!(rc, 2);
        assert!(err.contains("unrecognized arguments: --bogus"), "{err}");
    }
}
