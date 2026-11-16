//! `cut-haplotype START END FQ NAMES < SAM`: plants a deletion of the
//! reference bases START..END-1 (1-based) in the haplotype-1 reads of a
//! long-read slice (`tests/e2e/paralogs-paraphase-1.sh`).
//!
//! Each primary alignment tagged `HP:i:1` that overlaps the bases loses
//! them: a read across the gap is joined over it, one that starts or ends
//! inside keeps the part outside if that is 1 kb or more, else it is
//! dropped. Kept reads go to FQ as FASTQ, their RG, rq, HP and PS tags in the
//! comment (for `minimap2 -y`); every read that overlaps goes to NAMES, so
//! the caller can leave the originals out. Prints the count of joined reads,
//! and on stderr the three counts.

use super::{Fail, create, exactly, number};
use std::io::{BufRead, Write};
use std::path::Path;

/// What happened to one read.
#[derive(Debug, PartialEq, Eq)]
pub enum Cut {
    /// Not haplotype 1, not primary, or not over the bases.
    Skip,
    /// Joined across the gap, or trimmed: the FASTQ record.
    Spanning(String),
    Trimmed(String),
    Dropped,
}

/// The query offsets of the first base at or after `a` and `b`, and the
/// alignment's end (exclusive, 1-based).
fn offsets(
    pos: usize,
    cigar: &str,
    a: usize,
    b: usize,
) -> Result<(Option<usize>, Option<usize>, usize), String> {
    let (mut rf, mut qry) = (pos, 0usize);
    let (mut qa, mut qb) = (None, None);
    let mut n = 0usize;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + usize::from(c - b'0');
            continue;
        }
        match c {
            b'M' | b'=' | b'X' | b'D' | b'N' => {
                let consumes = matches!(c, b'M' | b'=' | b'X');
                for (x, q) in [(a, &mut qa), (b, &mut qb)] {
                    if q.is_none() && rf + n > x {
                        *q = Some(if consumes {
                            qry + x.saturating_sub(rf)
                        } else {
                            qry
                        });
                    }
                }
                rf += n;
                if consumes {
                    qry += n;
                }
            }
            b'I' | b'S' => qry += n,
            b'H' | b'P' => {}
            _ => return Err(format!("CIGAR '{cigar}' has an unknown operation")),
        }
        n = 0;
    }
    Ok((qa, qb, rf))
}

/// Cuts one SAM line.
pub fn cut(line: &str, a: usize, b: usize) -> Result<Cut, String> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 11 {
        return Err(format!("a SAM line with {} fields", f.len()));
    }
    let tags = &f[11..];
    let flag: u32 = f[1].parse().map_err(|_| format!("flag '{}'", f[1]))?;
    if !tags.contains(&"HP:i:1") || flag & 0x900 != 0 {
        return Ok(Cut::Skip);
    }
    let pos: usize = f[3].parse().map_err(|_| format!("POS '{}'", f[3]))?;
    let (qa, qb, end) = offsets(pos, f[5], a, b)?;
    if end <= a || pos >= b {
        return Ok(Cut::Skip);
    }
    let head = |s: &str| -> String {
        s.get(..qa.unwrap_or(0).min(s.len()))
            .unwrap_or("")
            .to_owned()
    };
    let tail = |s: &str| -> String {
        s.get(qb.unwrap_or(s.len()).min(s.len())..)
            .unwrap_or("")
            .to_owned()
    };
    let (left, lq) = if pos < a {
        (head(f[9]), head(f[10]))
    } else {
        (String::new(), String::new())
    };
    let (right, rq) = if end > b {
        (tail(f[9]), tail(f[10]))
    } else {
        (String::new(), String::new())
    };
    let keep: Vec<&str> = tags
        .iter()
        .copied()
        .filter(|t| {
            ["RG", "rq", "HP", "PS"]
                .iter()
                .any(|k| t.get(..2) == Some(k))
        })
        .collect();
    let record = format!(
        "@{}\t{}\n{left}{right}\n+\n{lq}{rq}\n",
        f[0],
        keep.join("\t")
    );
    Ok(if !left.is_empty() && !right.is_empty() {
        Cut::Spanning(record)
    } else if left.len() + right.len() >= 1000 {
        Cut::Trimmed(record)
    } else {
        Cut::Dropped
    })
}

pub fn run(
    args: &[String],
    input: &mut dyn BufRead,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<(), Fail> {
    let x = exactly(args, 4, "START END FQ NAMES")?;
    let (a, b) = (number(&x[0], "START")?, number(&x[1], "END")?);
    let (mut fq, mut names) = (create(Path::new(&x[2]))?, create(Path::new(&x[3]))?);
    let (mut spanning, mut trimmed, mut dropped) = (0u64, 0u64, 0u64);
    let mut line = String::new();
    while {
        line.clear();
        input.read_line(&mut line)? > 0
    } {
        let l = line.strip_suffix('\n').unwrap_or(&line);
        let c = cut(l, a, b).map_err(Fail::Run)?;
        if c != Cut::Skip {
            writeln!(names, "{}", l.split('\t').next().unwrap_or(""))?;
        }
        match c {
            Cut::Skip => {}
            Cut::Dropped => dropped += 1,
            Cut::Spanning(r) => {
                spanning += 1;
                fq.write_all(r.as_bytes())?;
            }
            Cut::Trimmed(r) => {
                trimmed += 1;
                fq.write_all(r.as_bytes())?;
            }
        }
    }
    fq.flush()?;
    names.flush()?;
    writeln!(
        err,
        "haplotype-1 reads joined across the deletion: {spanning}, trimmed: {trimmed}, dropped: {dropped}"
    )?;
    writeln!(out, "{spanning}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sam(name: &str, flag: u32, pos: usize, cigar: &str, seq: &str, tags: &str) -> String {
        format!(
            "{name}\t{flag}\tchr16\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{}\t{tags}",
            "#".repeat(seq.len())
        )
    }

    #[test]
    fn a_read_across_the_gap_is_joined_over_it() {
        // Bases 100..199 (1-based) of a read at 1 are reference 1..200; cut 101..150.
        let seq: String = (0..200).map(|i| char::from(b"ACGT"[i % 4])).collect();
        let line = sam("r", 0, 1, "200M", &seq, "RG:Z:x\tHP:i:1\tNM:i:0\tPS:i:5");
        let Cut::Spanning(r) = cut(&line, 101, 151).unwrap() else {
            panic!()
        };
        let want = format!(
            "@r\tRG:Z:x\tHP:i:1\tPS:i:5\n{}{}\n+\n{}\n",
            &seq[..100],
            &seq[150..],
            "#".repeat(150)
        );
        assert_eq!(r, want);
    }

    #[test]
    fn insertions_and_deletions_move_the_offsets() {
        // 10M 5I 10M 5D 10M from 1: ref 1-10, ins, 11-20, del 21-25, 26-35.
        assert_eq!(
            offsets(1, "10M5I10M5D10M", 13, 30).unwrap(),
            (Some(17), Some(29), 36)
        );
        // A cut point inside the deletion maps to the query base after it.
        assert_eq!(offsets(1, "10M5I10M5D10M", 22, 30).unwrap().0, Some(25));
        assert_eq!(offsets(1, "5S10M", 3, 8).unwrap(), (Some(7), Some(12), 11));
        assert!(offsets(1, "10Q", 3, 8).is_err());
    }

    #[test]
    fn others_are_skipped_short_ends_dropped_long_ends_trimmed() {
        let seq = "A".repeat(1200);
        let hp1 = "HP:i:1";
        assert_eq!(
            cut(&sam("r", 0, 1, "1200M", &seq, "HP:i:2"), 500, 600).unwrap(),
            Cut::Skip
        );
        assert_eq!(
            cut(&sam("r", 256, 1, "1200M", &seq, hp1), 500, 600).unwrap(),
            Cut::Skip
        );
        assert_eq!(
            cut(&sam("r", 2048, 1, "1200M", &seq, hp1), 500, 600).unwrap(),
            Cut::Skip
        );
        assert_eq!(
            cut(&sam("r", 0, 2000, "1200M", &seq, hp1), 500, 600).unwrap(),
            Cut::Skip
        );
        // Starts inside, keeps 1200 - 99 bases after 600: trimmed.
        assert!(matches!(
            cut(&sam("r", 0, 501, "1200M", &seq, hp1), 500, 600).unwrap(),
            Cut::Trimmed(_)
        ));
        // Starts inside, keeps 1200 - 399 bases after 900: dropped.
        assert_eq!(
            cut(&sam("r", 0, 501, "1200M", &seq, hp1), 500, 900).unwrap(),
            Cut::Dropped
        );
    }

    #[test]
    fn run_writes_the_fastq_the_names_and_the_count() {
        let d = crate::fixture::e2e::tests::scratch("cut");
        let seq = "C".repeat(3000);
        let input = [
            sam("span", 0, 1, "3000M", &seq, "HP:i:1"),
            sam("other", 0, 1, "3000M", &seq, "HP:i:2"),
            sam("drop", 0, 1500, "700M", &seq[..700], "HP:i:1"),
        ]
        .join("\n")
            + "\n";
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let p = |n: &str| d.join(n).to_str().unwrap().to_owned();
        let args = ["1100".to_owned(), "2100".to_owned(), p("fq"), p("names")];
        assert!(run(&args, &mut input.as_bytes(), &mut out, &mut err).is_ok());
        assert_eq!(String::from_utf8(out).unwrap(), "1\n");
        assert_eq!(
            std::fs::read_to_string(d.join("names")).unwrap(),
            "span\ndrop\n"
        );
        let fq = std::fs::read_to_string(d.join("fq")).unwrap();
        assert_eq!(fq.lines().nth(1).unwrap().len(), 2000);
        assert!(
            String::from_utf8(err)
                .unwrap()
                .contains("joined across the deletion: 1, trimmed: 0, dropped: 1")
        );
        std::fs::remove_dir_all(d).unwrap();
    }
}
