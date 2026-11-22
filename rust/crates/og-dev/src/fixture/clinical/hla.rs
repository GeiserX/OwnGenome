//! `hla-spike HLA_DAT OUT_DIR ALLELE...`: read pairs from the longest
//! genomic sequence in IMGT/HLA's `hla.dat` of each allele named, for the
//! spike-in control of `tests/e2e/hla-second-caller-1-bash.sh`. Pairs are
//! 150 bp ends of 350 bp fragments, one fragment every 5 bp, the strands
//! alternating: about 50x. They go to `OUT_DIR/R1.fq` and `OUT_DIR/R2.fq`.
//!
//! An entry counts for an allele when its `DE` name is the allele or starts
//! with the allele and a `:` (`HLA-B*57:01:01` takes `HLA-B*57:01:01:01`).
//! An allele with no entry, or none of 2500 bases or more, stops the kind.
//! Stdout: one line per allele with the entry taken, then `read pairs: N`.

use super::{Fail, bad, create};
use std::collections::HashMap;
use std::io::{BufRead, Write};

/// The reverse complement of `s`: `ACGTN` complemented, any other letter
/// kept, then reversed.
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

pub fn spike(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    let [dat, dir, wanted @ ..] = args else {
        return Err(Fail::Usage("takes HLA_DAT OUT_DIR ALLELE...".into()));
    };
    if wanted.is_empty() {
        return Err(Fail::Usage("takes at least one ALLELE".into()));
    }
    let mut best: HashMap<&str, (String, String)> = HashMap::new();
    let (mut name, mut seq, mut in_sq): (Option<String>, String, bool) =
        (None, String::new(), false);
    let input = og_io::open(dat).map_err(|e| bad(format!("{dat}: {e}")))?;
    for line in input.lines() {
        let line = line.map_err(|e| bad(format!("{dat}: {e}")))?;
        if let Some(rest) = line.strip_prefix("DE   ") {
            name = Some(rest.split(',').next().unwrap_or("").trim().to_string());
        } else if line.starts_with("SQ   ") {
            in_sq = true;
            seq.clear();
        } else if line.starts_with("//") {
            let s = seq.to_uppercase();
            if let Some(n) = name.as_deref().filter(|n| !n.is_empty()) {
                for w in wanted {
                    let prefix = format!("{w}:");
                    let longer = best
                        .get(w.as_str())
                        .is_none_or(|(_, b)| s.chars().count() > b.chars().count());
                    if (n == w || n.starts_with(&prefix)) && longer && !s.is_empty() {
                        best.insert(w, (n.to_string(), s.clone()));
                    }
                }
            }
            in_sq = false;
            name = None;
        } else if in_sq {
            seq.extend(line.chars().filter(|c| c.is_alphabetic()));
        }
    }
    let (mut r1, mut r2) = (
        create(&format!("{dir}/R1.fq"))?,
        create(&format!("{dir}/R2.fq"))?,
    );
    let qual = "I".repeat(150);
    let mut n = 0u64;
    for w in wanted {
        let Some((nm, s)) = best
            .get(w.as_str())
            .filter(|(_, s)| s.chars().count() >= 2500)
        else {
            r1.flush()?;
            r2.flush()?;
            return Err(bad(format!("no genomic sequence of {w} in {dat}")));
        };
        let chars: Vec<char> = s.chars().collect();
        writeln!(out, "{w}: {nm}, {} bp", chars.len())?;
        let mut start = 0;
        while start + 350 <= chars.len() {
            let frag = &chars[start..start + 350];
            let (mut a, mut b): (String, String) = (
                frag[..150].iter().collect(),
                revcomp(&frag[200..].iter().collect::<String>()),
            );
            if (start / 5) % 2 == 1 {
                std::mem::swap(&mut a, &mut b);
            }
            n += 1;
            write!(r1, "@spike{n}/1\n{a}\n+\n{qual}\n")?;
            write!(r2, "@spike{n}/2\n{b}\n+\n{qual}\n")?;
            start += 5;
        }
    }
    r1.flush()?;
    r2.flush()?;
    writeln!(out, "read pairs: {n}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{arg, dir, run};

    /// An `hla.dat`-like file: entries with DE, SQ and sequence lines.
    fn dat(entries: &[(&str, &str)]) -> String {
        let mut s = String::new();
        for (name, seq) in entries {
            s.push_str(&format!(
                "ID   X; SV 1; standard; DNA; HUM; {} BP.\n",
                seq.len()
            ));
            s.push_str(&format!("DE   {name}, Human MHC Class I sequence\n"));
            s.push_str("FT   exon 1..10\n");
            s.push_str(&format!("SQ   Sequence {} BP;\n", seq.len()));
            for chunk in seq.as_bytes().chunks(60) {
                s.push_str("     ");
                s.push_str(&String::from_utf8_lossy(chunk));
                s.push_str("        60\n");
            }
            s.push_str("//\n");
        }
        s
    }

    /// A made-up sequence of `n` bases, lower case as hla.dat writes them.
    fn bases(n: usize, seed: usize) -> String {
        (0..n)
            .map(|i| ['a', 'c', 'g', 't'][(i * 7 + seed + i / 13) % 4])
            .collect()
    }

    #[test]
    fn spike_takes_the_longest_matching_entry_and_tiles_it() {
        let d = dir("spike");
        let long = bases(2600, 1);
        std::fs::write(
            d.join("hla.dat"),
            dat(&[
                ("HLA-B*57:01:01:01", &bases(2550, 2)),
                ("HLA-B*57:01:01:02", &long),
                ("HLA-B*57:01:011", &bases(4000, 3)),
                ("HLA-B*07:02:01", &bases(3000, 4)),
            ]),
        )
        .expect("dat");
        let (code, out, err) = run(
            "hla-spike",
            &[
                &arg(&d, "hla.dat"),
                d.to_str().expect("utf8"),
                "HLA-B*57:01:01",
                "HLA-B*07:02:01",
            ],
        );
        assert_eq!(code, 0, "{err}");
        // (2600 - 350) / 5 + 1 = 451 and (3000 - 350) / 5 + 1 = 531 pairs.
        assert_eq!(
            out,
            "HLA-B*57:01:01: HLA-B*57:01:01:02, 2600 bp\nHLA-B*07:02:01: HLA-B*07:02:01, 3000 bp\nread pairs: 982\n"
        );
        let r1 = std::fs::read_to_string(d.join("R1.fq")).expect("r1");
        let r2 = std::fs::read_to_string(d.join("R2.fq")).expect("r2");
        let up = long.to_uppercase();
        let l1: Vec<&str> = r1.lines().collect();
        let l2: Vec<&str> = r2.lines().collect();
        assert_eq!(l1.len(), 982 * 4);
        assert_eq!(&l1[..4], ["@spike1/1", &up[..150], "+", &"I".repeat(150)]);
        assert_eq!(l2[1], super::revcomp(&up[200..350]));
        // The second fragment (start 5) has its strands swapped.
        assert_eq!(l1[5], super::revcomp(&up[205..355]));
        assert_eq!(l2[5], &up[5..155]);
    }

    #[test]
    fn spike_stops_on_an_allele_without_a_long_entry() {
        let d = dir("spike-short");
        std::fs::write(
            d.join("hla.dat"),
            dat(&[("HLA-A*01:01:01", &bases(2499, 5))]),
        )
        .expect("dat");
        let (code, _, err) = run(
            "hla-spike",
            &[
                &arg(&d, "hla.dat"),
                d.to_str().expect("utf8"),
                "HLA-A*01:01:01",
            ],
        );
        assert_eq!(code, 1);
        assert!(
            err.contains("no genomic sequence of HLA-A*01:01:01 in"),
            "{err}"
        );
    }
}
