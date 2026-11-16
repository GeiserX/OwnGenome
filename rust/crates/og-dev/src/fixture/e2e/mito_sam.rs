//! `mito-sam --fasta FA --header SAM --sample S --out OUT [--depth X]
//! [--deletion FIRST-LAST] [--level F]`: single-end 150-base chrM reads at
//! depth X (default 400) for step 20's case on large deletions
//! (`tests/e2e/mito-extended-1-deletions.sh`), the synth-sam of
//! `tests/test_mito_extended.py` with its seed (`random.Random(1)`).
//!
//! With `--deletion`, a share F (default 0.1) of the molecules lack the
//! bases FIRST..LAST (1-based). A read across the junction is written as an
//! aligner writes it: aligned on its longer side, the rest soft-clipped (a
//! forward read clipped on the right, a reverse read on the left); one with
//! a side under 30 bases is left out. The reads go under FA's `@HD` and
//! `@SQ` lines from SAM and one `@RG` of S. Prints how many reads it wrote
//! and how many are soft-clipped.

use super::{Fail, PyRandom, create, read};
use std::io::Write;
use std::path::Path;

const READ_LEN: usize = 150;

/// One read: 1-based position, CIGAR, flag and bases.
type Synth = (usize, String, u16, Vec<u8>);

/// The reads, sorted as Python sorts the tuples.
pub fn reads(seq: &[u8], depth: f64, deletion: Option<(usize, usize)>, level: f64) -> Vec<Synth> {
    let mut rnd = PyRandom::new(1);
    // int(len(seq) * depth / read_len): truncation of a float, as Python.
    let n = (seq.len() as f64 * depth / READ_LEN as f64) as usize;
    let mut out = Vec::with_capacity(n);
    let molecule = deletion.map(|(a, b)| ([&seq[..a - 1], &seq[b..]].concat(), a - 1));
    for _ in 0..n {
        if let (Some((a, b)), Some((mol, junction))) = (deletion, &molecule)
            && rnd.random() < level
        {
            let s = rnd.below(mol.len() - READ_LEN);
            let r = mol[s..s + READ_LEN].to_vec();
            if *junction <= s {
                out.push((
                    s + 1 + (b - a + 1),
                    format!("{READ_LEN}M"),
                    rnd.choice(&[0, 16]),
                    r,
                ));
                continue;
            }
            let left = junction - s;
            if left >= READ_LEN {
                out.push((s + 1, format!("{READ_LEN}M"), rnd.choice(&[0, 16]), r));
            } else if left.min(READ_LEN - left) < 30 {
                continue;
            } else if left >= READ_LEN - left {
                out.push((s + 1, format!("{left}M{}S", READ_LEN - left), 0, r));
            } else {
                out.push((b + 1, format!("{left}S{}M", READ_LEN - left), 16, r));
            }
            continue;
        }
        let s = rnd.below(seq.len() - READ_LEN);
        out.push((
            s + 1,
            format!("{READ_LEN}M"),
            rnd.choice(&[0, 16]),
            seq[s..s + READ_LEN].to_vec(),
        ));
    }
    out.sort();
    out
}

struct Opts {
    fasta: String,
    header: String,
    sample: String,
    out: String,
    depth: f64,
    deletion: Option<(usize, usize)>,
    level: f64,
}

fn parse(args: &[String]) -> Result<Opts, Fail> {
    let (mut fasta, mut header, mut sample, mut out) = (None, None, None, None);
    let (mut depth, mut deletion, mut level) = (400.0, None, 0.1);
    let float = |v: &str, k: &str| {
        v.parse::<f64>()
            .ok()
            .filter(|x| x.is_finite())
            .ok_or_else(|| Fail::Usage(format!("{k} '{v}' is not a number")))
    };
    let mut it = args.iter();
    while let Some(k) = it.next() {
        let v = it
            .next()
            .ok_or_else(|| Fail::Usage(format!("{k} takes a value")))?
            .clone();
        match k.as_str() {
            "--fasta" => fasta = Some(v),
            "--header" => header = Some(v),
            "--sample" => sample = Some(v),
            "--out" => out = Some(v),
            "--depth" => depth = float(&v, k)?,
            "--level" => level = float(&v, k)?,
            "--deletion" => {
                let span = v
                    .split_once('-')
                    .and_then(|(a, b)| Some((a.parse::<usize>().ok()?, b.parse::<usize>().ok()?)))
                    .filter(|&(a, b)| 1 < a && a <= b)
                    .ok_or_else(|| Fail::Usage(format!("--deletion '{v}' is not FIRST-LAST")))?;
                deletion = Some(span);
            }
            _ => return Err(Fail::Usage(format!("unknown option '{k}'"))),
        }
    }
    let need =
        |o: Option<String>, k: &str| o.ok_or_else(|| Fail::Usage(format!("{k} is required")));
    Ok(Opts {
        fasta: need(fasta, "--fasta")?,
        header: need(header, "--header")?,
        sample: need(sample, "--sample")?,
        out: need(out, "--out")?,
        depth,
        deletion,
        level,
    })
}

pub fn run(args: &[String], stdout: &mut dyn Write) -> Result<(), Fail> {
    let o = parse(args)?;
    let fasta = String::from_utf8_lossy(&read(Path::new(&o.fasta))?).into_owned();
    let seq: Vec<u8> = fasta
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if let Some((_, b)) = o.deletion
        && b + READ_LEN >= seq.len()
    {
        return Err(Fail::Run(format!(
            "--deletion ends past {} bases of chrM",
            seq.len()
        )));
    }
    if seq.len() <= READ_LEN {
        return Err(Fail::Run(format!(
            "{}: fewer than {READ_LEN} bases",
            o.fasta
        )));
    }
    let header = String::from_utf8_lossy(&read(Path::new(&o.header))?).into_owned();
    let reads = reads(&seq, o.depth, o.deletion, o.level);
    let mut w = create(Path::new(&o.out))?;
    for line in header.split_inclusive('\n') {
        if line.starts_with("@HD") || line.starts_with("@SQ") {
            w.write_all(line.as_bytes())?;
        }
    }
    let s = &o.sample;
    writeln!(w, "@RG\tID:{s}\tSM:{s}\tPL:ILLUMINA\tLB:{s}")?;
    for (i, (pos, cigar, flag, r)) in reads.iter().enumerate() {
        writeln!(
            w,
            "r{i}\t{flag}\tchrM\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t{}\tRG:Z:{s}",
            String::from_utf8_lossy(r),
            "I".repeat(r.len())
        )?;
    }
    w.flush()?;
    let clipped = reads.iter().filter(|r| r.1.contains('S')).count();
    writeln!(
        stdout,
        "{} reads, {clipped} soft-clipped at the junction",
        reads.len()
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture::e2e::tests::fnv;

    /// A made-up chrM: 16569 bases of `random.Random(5)`.
    fn chrm() -> Vec<u8> {
        PyRandom::new(5).bases(16569)
    }

    fn text(r: &[Synth]) -> String {
        r.iter()
            .map(|(p, c, f, s)| format!("{p} {c} {f} {}\n", String::from_utf8_lossy(s)))
            .collect()
    }

    // Pins of what tests/test_mito_extended.py's synth_reads gives on the
    // same chrM, with and without the common deletion.
    #[test]
    fn gives_the_pythons_reads() {
        let seq = chrm();
        let del = reads(&seq, 50.0, Some((8470, 13446)), 0.1);
        let wt = reads(&seq, 50.0, None, 0.1);
        assert_eq!((del.len(), wt.len()), PIN_COUNTS);
        assert_eq!(
            (fnv(text(&del).as_bytes()), fnv(text(&wt).as_bytes())),
            PIN_HASHES
        );
    }

    #[test]
    fn junction_reads_are_clipped_on_their_shorter_side() {
        let seq = chrm();
        let del = reads(&seq, 200.0, Some((8470, 13446)), 0.5);
        let clipped: Vec<&Synth> = del.iter().filter(|r| r.1.contains('S')).collect();
        assert!(clipped.len() > 10, "{}", clipped.len());
        for (pos, cigar, flag, _) in clipped {
            if cigar.ends_with('S') {
                assert_eq!(*flag, 0);
                assert!(*pos < 8470);
            } else {
                assert_eq!((*flag, *pos), (16, 13447));
            }
        }
        assert!(reads(&seq, 200.0, None, 0.5).iter().all(|r| r.1 == "150M"));
    }

    #[test]
    fn bad_options() {
        use crate::fixture::e2e::tests::run;
        assert_eq!(run("mito-sam", &["--fasta", "a"]).0, Some(2));
        assert_eq!(
            run(
                "mito-sam",
                &[
                    "--deletion",
                    "9-3",
                    "--fasta",
                    "a",
                    "--header",
                    "h",
                    "--sample",
                    "s",
                    "--out",
                    "o"
                ]
            )
            .0,
            Some(2)
        );
        assert_eq!(run("mito-sam", &["--bogus", "1"]).0, Some(2));
    }

    const PIN_COUNTS: (usize, usize) = (5520, 5523);
    const PIN_HASHES: (u64, u64) = (1_728_179_969_328_416_176, 8_620_160_239_877_651_625);
}
