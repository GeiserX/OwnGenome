//! `hba-tiles FASTA HEADER SAMPLE OUT`: forward 150-base reads every 5
//! bases over the alpha-globin window chr16:95001-235000 (FASTA, as
//! `samtools faidx REF chr16:95001-235000` writes it), under HEADER's lines,
//! for step 35's HBA case (`tests/e2e/paralog-checks-and-doc-fixes-1-hba.sh`).
//!
//! A read whose bases, or their reverse complement, occur more than once in
//! the window gets MAPQ 0, as an aligner gives a read that fits two places;
//! every other read MAPQ 60. A SAMPLE whose name ends in `3` loses every
//! other read that overlaps HBA2 (0-based 172875-173710), so HBA2 sits at
//! about half depth. Each read carries the first `@RG` of HEADER, if any.

use super::{Fail, create, exactly, read, revcomp};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

const START0: usize = 95_000;
const HBA2: (usize, usize) = (172_875, 173_710);
const L: usize = 150;

/// The SAM records (no header) of the window `seq`.
pub fn records(seq: &[u8], sample: &str, rg: Option<&str>) -> String {
    let mut seen: HashMap<&[u8], u32> = HashMap::new();
    for w in seq.windows(L) {
        *seen.entry(w).or_default() += 1;
    }
    let count = |r: &[u8]| seen.get(r).copied().unwrap_or(0);
    let half = sample.ends_with('3');
    let tag = rg
        .filter(|g| !g.is_empty())
        .map(|g| format!("\tRG:Z:{g}"))
        .unwrap_or_default();
    let q = "I".repeat(L);
    let mut out = String::new();
    let mut n = 0u64;
    for (i, r) in seq.windows(L).enumerate().step_by(5) {
        let s = START0 + i;
        if half && s < HBA2.1 && s + L > HBA2.0 {
            n += 1;
            if n % 2 == 1 {
                continue;
            }
        }
        let mapq = if count(r) + count(&revcomp(r)) > 1 {
            0
        } else {
            60
        };
        out += &format!(
            "syn{s}\t0\tchr16\t{}\t{mapq}\t{L}M\t*\t0\t0\t{}\t{q}{tag}\n",
            s + 1,
            String::from_utf8_lossy(r)
        );
    }
    out
}

pub fn run(args: &[String]) -> Result<(), Fail> {
    let a = exactly(args, 4, "FASTA HEADER SAMPLE OUT")?;
    let fasta = String::from_utf8_lossy(&read(Path::new(&a[0]))?).into_owned();
    let seq: Vec<u8> = fasta
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let header = read(Path::new(&a[1]))?;
    let text = String::from_utf8_lossy(&header);
    let rg = text
        .lines()
        .find(|l| l.starts_with("@RG"))
        .and_then(|l| l.split('\t').nth(1))
        .map(|f| f.get(3..).unwrap_or(""));
    let mut out = create(Path::new(&a[3]))?;
    out.write_all(&header)?;
    out.write_all(records(&seq, &a[2], rg).as_bytes())?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture::e2e::PyRandom;
    use crate::fixture::e2e::tests::fnv;

    /// 2000 random bases with a 300-base repeat, then its reverse complement.
    fn window() -> Vec<u8> {
        let mut s = PyRandom::new(9).bases(2000);
        let rep = s[100..400].to_vec();
        s[1000..1300].copy_from_slice(&rep);
        s[1500..1800].copy_from_slice(&revcomp(&rep));
        s
    }

    #[test]
    fn a_read_that_fits_twice_gets_mapq_0() {
        let recs = records(&window(), "S4", Some("rg1"));
        let mapq = |pos: usize| {
            recs.lines()
                .find(|l| l.split('\t').nth(3) == Some(&(START0 + pos + 1).to_string()))
                .map(|l| l.split('\t').nth(4).unwrap().to_owned())
        };
        assert_eq!(mapq(0).as_deref(), Some("60"));
        assert_eq!(mapq(150).as_deref(), Some("0")); // [150, 300) is inside the repeat
        assert_eq!(mapq(1000).as_deref(), Some("0"));
        assert_eq!(mapq(1550).as_deref(), Some("0")); // the reverse complement copy
        assert_eq!(mapq(600).as_deref(), Some("60"));
        assert!(recs.lines().all(|l| l.ends_with("\tRG:Z:rg1")));
        assert_eq!(recs.lines().count(), (2000 - L) / 5 + 1);
    }

    #[test]
    fn a_sample_ending_in_3_loses_every_other_read_over_hba2() {
        let seq = PyRandom::new(9).bases(80_000);
        let over = |r: &str| {
            r.lines()
                .filter(|l| {
                    let s: usize = l.split('\t').nth(3).unwrap().parse::<usize>().unwrap() - 1;
                    s < HBA2.1 && s + L > HBA2.0
                })
                .count()
        };
        let (four, three) = (
            records(&seq, "HG002hba4", None),
            records(&seq, "HG002hba3", None),
        );
        assert_eq!(over(&four), 196);
        assert_eq!(over(&three), 98);
        assert_eq!(four.lines().count() - three.lines().count(), 98);
        assert!(!four.contains("RG:Z"));
    }

    // A pin of what the inline Python of
    // tests/e2e/paralog-checks-and-doc-fixes-1-hba.sh writes on window().
    #[test]
    fn writes_the_pythons_records() {
        assert_eq!(
            fnv(records(&window(), "HG002hba3", Some("HG002")).as_bytes()),
            PIN
        );
    }

    const PIN: u64 = 12_081_885_455_233_452_406;
}
