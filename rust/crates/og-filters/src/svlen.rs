//! `og sv svlen-fix IN|-`: give AnnotSV the net length change of a
//! sequence-resolved SV record (replaces `bin/sv_svlen_fix.py`).
//!
//! Manta writes a deletion that also inserts bases with REF and ALT
//! sequences, `END` at the last deleted base and `SVLEN` as the deleted
//! length. AnnotSV takes the length from `SVLEN` and the coding overlap
//! from `POS-1..END`, so REF 137 bp with ALT 2 bp read as 137 deleted
//! coding bases: a frameshift and ACMG criterion 2E-1 for an in-frame
//! change. A record whose REF and one ALT are both plain sequence of
//! different lengths gets `SVLEN = len(ALT) - len(REF)` and, for a net
//! deletion, `END = POS - 1 + len(REF) - len(ALT)`; the old values stay in
//! `ORIG_SVLEN` and `ORIG_END`. Everything else passes through.
//!
//! A VCF with more than one sample column stops it (exit 2): step 05 reads
//! the first sample's genotype only.
//!
//! Fixes against the Python, each with a golden in `tests/parity`:
//! - an INFO field holding `SVLEN` or `END` twice: the Python saved the last
//!   value in `ORIG_*` but overwrote the first item, so the stale last one
//!   stayed beside the new value. The rewrite now replaces the item it read;
//! - a POS that is not a number on a record it rewrites: a traceback (exit
//!   1) before, now a message naming the line and exit 2.

use std::fmt;
use std::io::{self, BufRead, Write};

use crate::TextLines;

/// The two header lines the rewrite declares. The script name in them is
/// kept: the file they go into is a stable output.
const NEW_INFO: [&str; 2] = [
    "##INFO=<ID=ORIG_SVLEN,Number=1,Type=Integer,Description=\"SVLEN before sv_svlen_fix.py set it to the net length change len(ALT)-len(REF)\">",
    "##INFO=<ID=ORIG_END,Number=1,Type=Integer,Description=\"END before sv_svlen_fix.py set it to POS-1+len(REF)-len(ALT) for a net deletion, so AnnotSV reads the net length\">",
];

const USAGE: &str = "Usage: og sv svlen-fix IN.vcf[.gz]|- > OUT.vcf\n";

/// Why the input could not be fixed.
#[derive(Debug)]
pub enum SvlenError {
    Io(io::Error),
    /// More than one sample column: how many, and their names.
    MultiSample(Vec<String>),
    /// A record to rewrite whose POS is not a whole number: line and POS.
    BadPos(u64, String),
}

impl fmt::Display for SvlenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SvlenError::Io(e) => write!(f, "{e}"),
            SvlenError::MultiSample(names) => {
                let shown = names.iter().take(5).cloned().collect::<Vec<_>>().join(", ");
                let more = if names.len() > 5 { ", ..." } else { "" };
                write!(
                    f,
                    "the SV VCF holds {} samples ({shown}{more}); AnnotSV's input must hold one. Pick it first: bcftools view -s <sample>",
                    names.len()
                )
            }
            SvlenError::BadPos(line, pos) => write!(
                f,
                "{line}: POS {} is not a whole number",
                crate::py_repr(pos)
            ),
        }
    }
}

impl From<io::Error> for SvlenError {
    fn from(e: io::Error) -> Self {
        SvlenError::Io(e)
    }
}

/// Plain sequence: `^[ACGTNacgtn]+$`.
fn is_seq(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|c| b"ACGTNacgtn".contains(&c))
}

/// The net length change of a record whose REF and one ALT are plain
/// sequence of different lengths; `None` for every other record.
pub fn net_change(ref_: &str, alt: &str) -> Option<i64> {
    if alt.contains(',') || !is_seq(ref_) || !is_seq(alt) || ref_.len() == alt.len() {
        return None;
    }
    i64::try_from(alt.len())
        .ok()
        .zip(i64::try_from(ref_.len()).ok())
        .map(|(a, r)| a - r)
}

/// The rewritten INFO of a record, or `None` when it stays as it is.
/// `fields` are the record's columns (8 or more).
fn fix(fields: &[&str], line: u64) -> Result<Option<String>, SvlenError> {
    let Some(net) = net_change(fields[3], fields[4]) else {
        return Ok(None);
    };
    let pos = pycompat::num::parse_int(fields[1])
        .map_err(|_| SvlenError::BadPos(line, fields[1].to_owned()))?;
    let info = fields[7];
    let mut items: Vec<String> = if info.is_empty() || info == "." {
        Vec::new()
    } else {
        info.split(';').map(str::to_owned).collect()
    };
    // The value of a key is its last `key=value` item, as a dict of the
    // items reads it; `at` is that item, the one a rewrite replaces.
    let lookup = |items: &[String], key: &str| -> Option<(usize, String)> {
        items.iter().enumerate().rev().find_map(|(i, it)| {
            let (k, v) = it.split_once('=')?;
            (k == key).then(|| (i, v.to_owned()))
        })
    };
    let net_s = net.to_string();
    let svtype =
        lookup(&items, "SVTYPE").map_or(if net < 0 { "DEL" } else { "INS" }.to_owned(), |(_, v)| v);
    let mut changed: Vec<(&str, String)> = Vec::new();
    if lookup(&items, "SVLEN").map(|(_, v)| v).as_deref() != Some(net_s.as_str()) {
        changed.push(("SVLEN", net_s));
    }
    if net < 0 && svtype == "DEL" {
        let ref_len = i128::from(fields[3].len() as u64);
        let alt_len = i128::from(fields[4].len() as u64);
        let end = (pos - 1 + ref_len - alt_len).to_string();
        if lookup(&items, "END").map(|(_, v)| v).as_deref() != Some(end.as_str()) {
            changed.push(("END", end));
        }
    }
    if changed.is_empty() {
        return Ok(None);
    }
    for (key, value) in changed {
        match lookup(&items, key) {
            Some((at, old)) => {
                items.push(format!("ORIG_{key}={old}"));
                items[at] = format!("{key}={value}");
            }
            None => items.push(format!("{key}={value}")),
        }
    }
    Ok(Some(items.join(";")))
}

/// Copies the VCF from `input` to `out`, fixed; the number of records it
/// rewrote.
pub fn run<R: BufRead>(mut lines: TextLines<R>, mut out: impl Write) -> Result<usize, SvlenError> {
    let mut line = String::new();
    let mut seen = [false; 2];
    let mut n = 0;
    let mut number = 0u64;
    while lines.read_line(&mut line)? {
        number += 1;
        let body = line.strip_suffix('\n').unwrap_or(&line);
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                let samples: Vec<String> = body.split('\t').skip(9).map(str::to_owned).collect();
                if samples.len() > 1 {
                    return Err(SvlenError::MultiSample(samples));
                }
                for (h, was) in NEW_INFO.iter().zip(seen) {
                    if !was {
                        writeln!(out, "{h}")?;
                    }
                }
            }
            for (h, was) in NEW_INFO.iter().zip(seen.iter_mut()) {
                *was |= body == *h;
            }
            out.write_all(line.as_bytes())?;
            continue;
        }
        let fields: Vec<&str> = body.split('\t').collect();
        let new_info = if fields.len() >= 8 {
            fix(&fields, number)?
        } else {
            None
        };
        match new_info {
            Some(info) => {
                n += 1;
                let mut rest = fields;
                rest[7] = &info;
                out.write_all(rest.join("\t").as_bytes())?;
            }
            None => out.write_all(body.as_bytes())?,
        }
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(n)
}

/// `og sv svlen-fix IN|-`: exit 0, or 2 for a usage error or an input it
/// cannot fix.
pub fn cli(args: &[String]) -> u8 {
    let [path] = args else {
        eprint!("{USAGE}");
        return 2;
    };
    if path == "-h" || path == "--help" {
        print!("{USAGE}");
        return 0;
    }
    let opened = if path == "-" {
        og_io::open(path).map(TextLines::stdin)
    } else {
        og_io::open(path).map(TextLines::file)
    };
    let input = match opened {
        Ok(i) => i,
        Err(e) => {
            eprintln!("ERROR: {}", crate::os_error(&e, path));
            return 2;
        }
    };
    let stdout = io::stdout().lock();
    match run(input, io::BufWriter::with_capacity(1 << 16, stdout)) {
        Ok(n) => {
            eprintln!(
                "SV length check: {n} sequence-resolved records set to their net length change"
            );
            0
        }
        Err(SvlenError::MultiSample(names)) => {
            eprintln!("ERROR: {}", SvlenError::MultiSample(names));
            2
        }
        Err(e) => {
            eprintln!("ERROR: {path}:{e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixed(vcf: &str) -> (String, usize) {
        let mut out = Vec::new();
        let n = run(TextLines::file(vcf.as_bytes()), &mut out).unwrap();
        (String::from_utf8(out).unwrap(), n)
    }

    const HEAD: &str =
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

    #[test]
    fn net_change_only_for_plain_sequence_of_two_lengths() {
        assert_eq!(net_change("ACGT", "A"), Some(-3));
        assert_eq!(net_change("A", "ACG"), Some(2));
        assert_eq!(net_change("A", "C"), None);
        assert_eq!(net_change("A", "<DEL>"), None);
        assert_eq!(net_change("A", "AC,AG"), None);
        assert_eq!(net_change("A", "N]chr1:5]"), None);
        assert_eq!(net_change("", "A"), None);
    }

    #[test]
    fn a_deletion_that_inserts_gets_its_net_change() {
        let r = "A".repeat(137);
        let vcf = format!(
            "{HEAD}chr20\t1005400\tm\t{r}\tAC\t.\tPASS\tEND=1005536;SVTYPE=DEL;SVLEN=-136\tGT\t0/1\n"
        );
        let (out, n) = fixed(&vcf);
        assert_eq!(n, 1);
        assert!(
            out.contains("\tEND=1005534;SVTYPE=DEL;SVLEN=-135;ORIG_SVLEN=-136;ORIG_END=1005536\t"),
            "{out}"
        );
        assert!(out.contains(NEW_INFO[0]) && out.contains(NEW_INFO[1]));
        // A second run changes nothing and declares nothing twice.
        let (again, n) = fixed(&out);
        assert_eq!((again.as_str(), n), (out.as_str(), 0));
    }

    #[test]
    fn red_first_the_end_is_one_base_short_of_annotsvs_default() {
        // AnnotSV's own END is POS + len(REF) - len(ALT); ours is one less.
        let vcf = format!("{HEAD}c\t100\t.\tACGTA\tA\t.\tPASS\t.\n");
        let (out, _) = fixed(&vcf);
        assert!(out.ends_with("\tSVLEN=-4;END=103\n"), "{out}");
    }

    #[test]
    fn a_duplicated_key_rewrites_the_item_it_read() {
        let vcf = format!("{HEAD}c\t100\t.\tACGTA\tA\t.\tPASS\tSVLEN=-1;SVLEN=-2;END=103\n");
        let (out, _) = fixed(&vcf);
        assert!(
            out.ends_with("\tSVLEN=-1;SVLEN=-4;END=103;ORIG_SVLEN=-2\n"),
            "{out}"
        );
    }

    #[test]
    fn short_lines_and_symbolic_records_pass_through() {
        let vcf = format!("{HEAD}c\t1\n\nc\t5\t.\tA\t<DEL>\t.\tPASS\tSVLEN=-9\n");
        let (out, n) = fixed(&vcf);
        assert_eq!(n, 0);
        assert!(out.ends_with("c\t1\n\nc\t5\t.\tA\t<DEL>\t.\tPASS\tSVLEN=-9\n"));
    }

    #[test]
    fn several_samples_and_a_bad_pos_are_errors() {
        let two = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n";
        let e = run(TextLines::file(two.as_bytes()), Vec::new()).unwrap_err();
        assert_eq!(
            e.to_string(),
            "the SV VCF holds 2 samples (a, b); AnnotSV's input must hold one. Pick it first: bcftools view -s <sample>"
        );
        let bad = format!("{HEAD}c\tx\t.\tACG\tA\t.\t.\t.\n");
        let e = run(TextLines::file(bad.as_bytes()), Vec::new()).unwrap_err();
        assert_eq!(e.to_string(), "3: POS 'x' is not a whole number");
    }
}
