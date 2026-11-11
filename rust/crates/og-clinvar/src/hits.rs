//! One row per hit of step 06's hits VCF (`og clinvar hits`).
//!
//! Columns, tab-separated, no header: STARS CHROM POS REF ALT GENOTYPE GENE
//! SIGNIFICANCE REVIEW_STATUS COMBINED_ALLELE DP GQ LOW_CONFIDENCE.
//!
//! - STARS from [`ReviewStatus`]; GENOTYPE from [`zygosity`]; GENE the
//!   GENEINFO symbols, comma-joined (`.` when none); SIGNIFICANCE and
//!   REVIEW_STATUS are CLNSIG and CLNREVSTAT with `_` as spaces (`.` when
//!   absent); COMBINED_ALLELE from [`combined_note`].
//! - DP and GQ are FORMAT/DP and FORMAT/GQ when they are whole numbers,
//!   else `.`; LOW_CONFIDENCE is `yes` when either is below the call floor
//!   (`##PGP_CALL_FLOOR=GQ=..,DP=..`, GQ 30 and DP 10 without the line).

use std::io::{self, BufRead, Write};

use crate::text::{is_digits, nth, split};
use crate::{CallFloor, ReviewStatus, combined_note, zygosity};

/// `_` as spaces, or `.` for an empty value.
fn words(v: &[u8]) -> Vec<u8> {
    if v.is_empty() {
        return b".".to_vec();
    }
    v.iter()
        .map(|&c| if c == b'_' { b' ' } else { c })
        .collect()
}

/// The row of one hit record, under `floor`.
pub fn row(line: &[u8], floor: &CallFloor) -> Vec<u8> {
    let f = split(line, b'\t');
    let (mut geneinfo, mut clnsig, mut rev): (&[u8], &[u8], &[u8]) = (b"", b"", b"");
    let (mut joined, mut part_of, mut one): (&[u8], &[u8], &[u8]) = (b"", b"", b"");
    for kv in split(nth(&f, 8), b';') {
        let Some(eq) = kv.iter().position(|&c| c == b'=') else {
            continue;
        };
        let v = &kv[eq + 1..];
        match &kv[..eq] {
            b"GENEINFO" => geneinfo = v,
            b"CLNSIG" => clnsig = v,
            b"CLNREVSTAT" => rev = v,
            b"PGP_JOINED_FROM" => joined = v,
            b"PGP_PART_OF" => part_of = v,
            b"PGP_ONE_ALLELE" => one = v,
            _ => {}
        }
    }
    let stars = ReviewStatus::parse(rev).stars();
    let genes: Vec<&[u8]> = split(geneinfo, b'|')
        .into_iter()
        .filter(|g| !g.is_empty())
        .map(|g| nth(&split(g, b':'), 1))
        .collect();
    let gene = genes.join(&b',');
    let gene = if gene.is_empty() { b".".to_vec() } else { gene };
    let sample = split(nth(&f, 10), b':');
    let gt = nth(&sample, 1);
    let (mut dp, mut gq): (&[u8], &[u8]) = (b".", b".");
    for (i, k) in split(nth(&f, 9), b':').into_iter().enumerate() {
        let v = nth(&sample, i + 1);
        if k == b"DP" && is_digits(v) {
            dp = v;
        } else if k == b"GQ" && is_digits(v) {
            gq = v;
        }
    }
    let low = floor.low((gq != b".").then_some(gq), (dp != b".").then_some(dp));
    [
        stars.to_string().as_bytes(),
        nth(&f, 1),
        nth(&f, 2),
        nth(&f, 4),
        nth(&f, 5),
        &zygosity(gt),
        &gene,
        &words(clnsig),
        &words(rev),
        &combined_note(joined, part_of, one, gt),
        dp,
        gq,
        if low { b"yes" } else { b"no" },
    ]
    .join(&b'\t')
}

/// The rows of a hits VCF. A `##PGP_CALL_FLOOR=` line sets the floor for
/// the records after it.
pub fn hits(mut input: impl BufRead, out: &mut impl Write) -> io::Result<()> {
    let mut floor = CallFloor::default();
    let mut line = Vec::new();
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Ok(());
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.starts_with(b"##PGP_CALL_FLOOR=") {
            floor.read_line(&line);
        } else if !line.starts_with(b"#") {
            out.write_all(&row(&line, &floor))?;
            out.write_all(b"\n")?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(line: &str) -> String {
        String::from_utf8(row(line.as_bytes(), &CallFloor::default())).unwrap()
    }

    #[test]
    fn a_hit_in_thirteen_columns() {
        assert_eq!(
            r(
                "chr1\t10\t5\tC\tT\t.\tPASS\tGENEINFO=A:1|B:2;CLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT:GQ:DP\t0/1:20:30"
            ),
            "3\tchr1\t10\tC\tT\thet\tA,B\tLikely pathogenic\treviewed by expert panel\t.\t30\t20\tyes"
        );
    }

    #[test]
    fn missing_values_and_an_empty_gene_segment() {
        assert_eq!(
            r("chr1\t10\t.\tC\tT\t.\tPASS\t.\tGT\t1/."),
            "0\tchr1\t10\tC\tT\tone copy, other copy not called ref (1/.)\t.\t.\t.\t.\t.\t.\tno"
        );
        assert_eq!(
            r("c\t1\t.\tC\tT\t.\t.\tGENEINFO=A:1||:3\tGT:DP\t1/1:x"),
            "0\tc\t1\tC\tT\thom\tA,\t.\t.\t.\t.\t.\tno"
        );
    }
}
