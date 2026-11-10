//! A FASTA with every base inside a BED row set to `N`, everything else
//! copied byte for byte: headers, line lengths, line endings, lower case.
//!
//! One pass, one line in memory: each record's rows are merged into
//! disjoint ascending intervals, and a cursor walks them as the lines go by,
//! so the cost is lines plus rows, never lines times rows.

use crate::bed::{BedError, Rows};
use std::collections::HashMap;
use std::fmt;
use std::io::{self, BufRead, Write};

/// What [`mask_fasta`] did, printed as one line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Summary {
    /// Rows that covered at least one base of a sequence in the FASTA.
    pub rows_masked: usize,
    /// Every row of the BED.
    pub rows_total: usize,
    /// Bases inside a row, each counted once however many rows cover it
    /// (and once per record, when a name appears twice).
    pub bases: u64,
    /// Sequences of the BED that no FASTA header names, sorted.
    pub missing: Vec<String>,
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let missing = if self.missing.is_empty() {
            "none".to_owned()
        } else {
            self.missing.join(", ")
        };
        write!(
            f,
            "masked {} of {} rows ({} bases); rows on sequences the FASTA lacks: {missing}",
            self.rows_masked, self.rows_total, self.bases
        )
    }
}

/// Why the FASTA could not be masked.
#[derive(Debug)]
pub enum MaskError {
    Io(io::Error),
    /// The FASTA cannot be read as lines of bases: its line (1-based) and why.
    Fasta {
        line: u64,
        why: &'static str,
    },
}

impl fmt::Display for MaskError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MaskError::Io(e) => write!(f, "{e}"),
            MaskError::Fasta { line, why } => write!(f, "line {line}: {why}"),
        }
    }
}

impl std::error::Error for MaskError {}

impl From<io::Error> for MaskError {
    fn from(e: io::Error) -> Self {
        MaskError::Io(e)
    }
}

/// `bytes.split()[0]` of a header after its `>`: the name ends at ASCII
/// whitespace.
fn header_name(after_gt: &[u8]) -> Option<&[u8]> {
    after_gt
        .split(|c| matches!(c, b' ' | b'\t' | b'\n' | b'\r' | b'\x0b' | b'\x0c'))
        .find(|w| !w.is_empty())
}

/// The record being copied: its intervals, the cursor into them, and the
/// position of the next line's first base.
struct Record<'a> {
    name: Option<&'a str>,
    iv: &'a [(u64, u64)],
    next: usize,
    pos: u64,
}

/// Records how long `rec` was, for the rows-masked count.
fn close<'a>(rec: &Record<'a>, seen: &mut HashMap<&'a str, u64>) {
    if let Some(name) = rec.name {
        let len = seen.entry(name).or_insert(0);
        *len = (*len).max(rec.pos);
    }
}

/// Copies `src` to `out` with every base inside a row of `rows` set to `N`.
///
/// A line starting with `>` is a header; the first word after the `>` names
/// the record, and the record's bases are the bytes of its following lines
/// without their line ending (`\n` or `\r\n`). A carriage return anywhere
/// else (old Mac line endings) is refused: it would be counted as a base and
/// shift every coordinate after it. Lines before the first header are copied
/// as they are. A name that appears twice is masked in both records.
pub fn mask_fasta(
    src: impl BufRead,
    rows: &Rows,
    mut out: impl Write,
) -> Result<Summary, MaskError> {
    let mut src = src;
    let mut line: Vec<u8> = Vec::with_capacity(256);
    let mut rec = Record {
        name: None,
        iv: &[],
        next: 0,
        pos: 0,
    };
    // The longest record seen for each name that has rows.
    let mut seen: HashMap<&str, u64> = HashMap::new();
    let mut bases: u64 = 0;
    let mut n: u64 = 0;
    loop {
        line.clear();
        if src.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        n += 1;
        let body_len = body_len(&line).ok_or(MaskError::Fasta {
            line: n,
            why: "a carriage return that does not end the line (old Mac line endings?); convert the file to \\n or \\r\\n line endings",
        })?;
        if line.first() == Some(&b'>') {
            close(&rec, &mut seen);
            let name = header_name(&line[1..]).ok_or(MaskError::Fasta {
                line: n,
                why: "a header without a sequence name",
            })?;
            // The name as the key `rows` holds, so it lives as long as `rows`.
            let name = std::str::from_utf8(name).ok().and_then(|s| rows.key(s));
            rec = Record {
                name,
                iv: name.map_or(&[], |s| rows.merged(s)),
                next: 0,
                pos: 0,
            };
            out.write_all(&line)?;
            continue;
        }
        let start = rec.pos;
        let end = start + body_len as u64;
        while rec.next < rec.iv.len() && rec.iv[rec.next].1 <= start {
            rec.next += 1;
        }
        let mut k = rec.next;
        while k < rec.iv.len() && rec.iv[k].0 < end {
            let (s, e) = rec.iv[k];
            let (a, b) = (s.max(start), (e + 1).min(end));
            // Both are inside the line, so they fit in usize.
            line[(a - start) as usize..(b - start) as usize].fill(b'N');
            bases += b - a;
            if e > end {
                break;
            }
            k += 1;
        }
        rec.pos = end;
        out.write_all(&line)?;
    }
    close(&rec, &mut seen);
    out.flush()?;

    let mut rows_masked = 0;
    let mut missing = Vec::new();
    for name in rows.names() {
        match seen.get(name) {
            Some(&len) => {
                let seq = rows.seq(name).expect("a name of rows");
                rows_masked += seq.rows.iter().filter(|&&(s, e)| s < e && s < len).count();
            }
            None => missing.push(name.to_owned()),
        }
    }
    missing.sort_unstable();
    Ok(Summary {
        rows_masked,
        rows_total: rows.total(),
        bases,
        missing,
    })
}

/// The length of a line without its ending (`\n`, `\r\n` or none at the
/// end of the file), or `None` when it holds any other carriage return.
fn body_len(line: &[u8]) -> Option<usize> {
    let body = line
        .strip_suffix(b"\r\n")
        .or_else(|| line.strip_suffix(b"\n"))
        .unwrap_or(line);
    (!body.contains(&b'\r')).then_some(body.len())
}

/// `og setup mask-fasta SRC BED DEST`: writes DEST (atomically: it appears
/// only once whole) and prints the [`Summary`]. Exit 0, 1 on an unreadable
/// or refused input, 2 on a wrong command line.
pub fn cli(args: &[String]) -> u8 {
    let [src, bed, dest] = args else {
        eprintln!("usage: og setup mask-fasta SRC BED DEST");
        return 2;
    };
    let fail = |what: &str, e: &dyn fmt::Display| {
        eprintln!("og setup mask-fasta: {what}: {e}");
        1
    };
    let rows = match og_io::open(bed).map_err(BedError::Io).and_then(Rows::read) {
        Ok(r) => r,
        Err(e) => return fail(bed, &e),
    };
    let input = match std::fs::File::open(src) {
        Ok(f) => io::BufReader::with_capacity(1 << 20, f),
        Err(e) => return fail(src, &e),
    };
    let mut out = match og_io::AtomicFile::create(dest) {
        Ok(o) => o,
        Err(e) => return fail(dest, &e),
    };
    let summary = match mask_fasta(input, &rows, &mut out) {
        Ok(s) => s,
        Err(MaskError::Io(e)) => return fail(&format!("{src} -> {dest}"), &e),
        Err(e) => return fail(src, &e),
    };
    if let Err(e) = out.commit() {
        return fail(dest, &e);
    }
    println!("{summary}");
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(fa: &[u8], bed: &str) -> Result<(Vec<u8>, String), MaskError> {
        let rows = Rows::read(bed.as_bytes()).expect("bed");
        let mut out = Vec::new();
        let s = mask_fasta(fa, &rows, &mut out)?;
        Ok((out, s.to_string()))
    }

    #[test]
    fn the_fixitfelix_example() {
        // s1 is ACGTacgtACGTACGT: 0-based [5, 9) and [15, 16) are N.
        let (out, s) = run(
            b">s1 first\nACGTacg\ntACGTAC\nGT\n>s2\nAAAAAAA\nCC\n",
            "browser position s1:1-3\n#seq\tstart\tend\ns1\t5\t9\ns1\t15\t16\nmissing\t0\t10\n",
        )
        .unwrap();
        assert_eq!(out, b">s1 first\nACGTaNN\nNNCGTAC\nGN\n>s2\nAAAAAAA\nCC\n");
        assert_eq!(
            s,
            "masked 2 of 3 rows (5 bases); rows on sequences the FASTA lacks: missing"
        );
    }

    #[test]
    fn crlf_keeps_its_line_endings() {
        let (out, s) = run(b">s1\r\nACGT\r\nACGT\r\n", "s1 3 5\n").unwrap();
        assert_eq!(out, b">s1\r\nACGN\r\nNCGT\r\n");
        assert_eq!(
            s,
            "masked 1 of 1 rows (2 bases); rows on sequences the FASTA lacks: none"
        );
    }

    #[test]
    fn last_line_without_a_newline() {
        let (out, _) = run(b">s1\nACGT\nAC", "s1 0 6\n").unwrap();
        assert_eq!(out, b">s1\nNNNN\nNN");
    }

    #[test]
    fn lone_carriage_returns_are_refused() {
        for fa in [
            &b">s1\rACGT\rACGT\r"[..],
            b">s1\nAC\rGT\n",
            b">s1\nACGT\r",
            b">s1\nACGT\r\r\n",
        ] {
            let e = run(fa, "s1 0 1\n").unwrap_err();
            assert!(e.to_string().contains("carriage return"), "{fa:?}: {e}");
        }
    }

    #[test]
    fn overlapping_rows_count_each_base_once() {
        let (out, s) = run(b">s1\nAAAAAAAAAA\n", "s1 1 4\ns1 2 6\ns1 2 3\n").unwrap();
        assert_eq!(out, b">s1\nANNNNNAAAA\n");
        assert_eq!(
            s,
            "masked 3 of 3 rows (5 bases); rows on sequences the FASTA lacks: none"
        );
    }

    #[test]
    fn rows_past_the_end_and_empty_rows_mask_nothing() {
        let (out, s) = run(b">s1\nAAAA\nAA\n", "s1 4 100\ns1 6 9\ns1 2 2\n").unwrap();
        assert_eq!(out, b">s1\nAAAA\nNN\n");
        assert_eq!(
            s,
            "masked 1 of 3 rows (2 bases); rows on sequences the FASTA lacks: none"
        );
    }

    #[test]
    fn a_row_spanning_many_lines_and_a_repeated_name() {
        let (out, s) = run(
            b"x\n>s1 a\nAAA\nAAA\nAAA\n>s2\nCC\n>s1 b\nAAA\n",
            "s1 2 7\n",
        )
        .unwrap();
        assert_eq!(out, b"x\n>s1 a\nAAN\nNNN\nNAA\n>s2\nCC\n>s1 b\nAAN\n");
        assert_eq!(
            s,
            "masked 1 of 1 rows (6 bases); rows on sequences the FASTA lacks: none"
        );
    }

    #[test]
    fn empty_bed_copies_the_file() {
        let fa = b">s1\nacgtNNNN\n>s2\n\n";
        let (out, s) = run(fa, "").unwrap();
        assert_eq!(out, fa);
        assert_eq!(
            s,
            "masked 0 of 0 rows (0 bases); rows on sequences the FASTA lacks: none"
        );
    }

    #[test]
    fn a_header_without_a_name_is_refused() {
        let e = run(b">s1\nA\n> \nC\n", "s1 0 1\n").unwrap_err();
        assert_eq!(e.to_string(), "line 3: a header without a sequence name");
    }

    #[test]
    fn missing_sequences_are_sorted() {
        let (_, s) = run(b">s1\nA\n", "zz 0 1\nab 0 1\ns1 0 1\n").unwrap();
        assert_eq!(
            s,
            "masked 1 of 3 rows (1 bases); rows on sequences the FASTA lacks: ab, zz"
        );
    }
}
