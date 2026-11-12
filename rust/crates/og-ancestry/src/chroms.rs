//! `chroms`: the chromosomes to phase and paint, those with N genotyped
//! SNVs or more (default 100) spanning X cM or more (default 5). RFMix
//! stops with a segmentation fault on a span of a few cM; on a whole
//! genome every autosome spans more than 40.

use crate::{Args, Fail, Result, Spec};
use std::io::Write;
use std::path::Path;

pub(crate) const SPEC: Spec = Spec {
    verb: "chroms",
    usage: "usage: og ancestry-fine chroms --counts counts.tsv --maps DIR [--min-sites N] [--min-cm X]\n",
    options: &["--counts", "--maps", "--min-sites", "--min-cm"],
    required: &["--counts", "--maps"],
};

/// A PLINK `.map` (chromosome, id, cM, position) as (position, cM),
/// sorted.
pub fn read_map(path: &Path) -> Result<Vec<(i64, f64)>> {
    let mut out = Vec::new();
    let mut n = 0u64;
    crate::for_each_line(path, |line| {
        n += 1;
        let c: Vec<&str> = crate::split_ws(line).collect();
        if c.len() >= 4 {
            let bad = |what: &str, v: &str| {
                Fail(format!(
                    "{}: line {n}: {what}: {}",
                    path.display(),
                    crate::repr(v)
                ))
            };
            let pos =
                crate::int(c[3]).ok_or_else(|| bad("position is not a whole number", c[3]))?;
            let cm =
                pycompat::num::parse_float(c[2]).map_err(|_| bad("cM is not a number", c[2]))?;
            out.push((pos, cm));
        }
        Ok(())
    })?;
    out.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
    Ok(out)
}

/// The genetic position of `pos`, linear between the map's markers.
pub fn cm_at(map: &[(i64, f64)], pos: i64) -> f64 {
    let (Some(&(p_first, c_first)), Some(&(p_last, c_last))) = (map.first(), map.last()) else {
        return 0.0;
    };
    if pos <= p_first {
        return c_first;
    }
    if pos >= p_last {
        return c_last;
    }
    // The last marker at or before pos, and the one after it.
    let hi = map.partition_point(|&(p, _)| p <= pos);
    let ((p0, c0), (p1, c1)) = (map[hi - 1], map[hi]);
    if p1 > p0 {
        c0 + (c1 - c0) * (pos - p0) as f64 / (p1 - p0) as f64
    } else {
        c0
    }
}

/// One row of `counts.tsv`: chromosome, genotyped SNVs, first and last
/// position.
pub type Count = (String, i64, i64, i64);

pub fn read_counts(path: &Path) -> Result<Vec<Count>> {
    let mut out = Vec::new();
    let mut n = 0u64;
    crate::for_each_line(path, |line| {
        n += 1;
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() >= 4 && c[0].starts_with("chr") {
            let num = |v: &str| {
                crate::int(v).ok_or_else(|| {
                    Fail(format!(
                        "{}: line {n}: not a whole number: {}",
                        path.display(),
                        crate::repr(v)
                    ))
                })
            };
            out.push((c[0].to_owned(), num(c[1])?, num(c[2])?, num(c[3])?));
        }
        Ok(())
    })?;
    Ok(out)
}

/// The chromosomes of `counts` to phase and paint, in chromosome order.
pub fn paint_chroms(
    counts: &[Count],
    maps: &Path,
    min_sites: i64,
    min_cm: f64,
) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for (chrom, n, first, last) in counts {
        let path = maps.join(format!("{chrom}.map"));
        if *n < min_sites || !path.is_file() {
            continue;
        }
        let map = read_map(&path)?;
        if cm_at(&map, *last) - cm_at(&map, *first) >= min_cm {
            out.push(chrom.clone());
        }
    }
    out.sort_by(|a, b| crate::chrom_key(a).cmp(&crate::chrom_key(b)));
    Ok(out)
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let counts = read_counts(Path::new(a.req("--counts")))?;
    let min_sites = a.int("--min-sites", 100).unwrap_or(100);
    let min_cm = a.float("--min-cm", 5.0).unwrap_or(5.0);
    let chroms = paint_chroms(&counts, Path::new(a.req("--maps")), min_sites, min_cm)?;
    writeln!(out, "{}", chroms.join(" ")).map_err(|e| Fail(format!("cannot write stdout: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    #[test]
    fn a_chromosome_needs_its_snvs_and_its_span() {
        let d = Dir::new("chroms");
        for c in ["chr20", "chr22", "chr19"] {
            d.put(
                &format!("maps/{c}.map"),
                &format!("{c} . 0.0 1000000\n{c} . 10.0 11000000\n"),
            );
        }
        let counts: Vec<Count> = vec![
            ("chr20".into(), 500, 2_000_000, 8_000_000),
            ("chr22".into(), 500, 2_000_000, 4_000_000),
            ("chr21".into(), 999, 1, 40_000_000),
            ("chr19".into(), 50, 1_000_000, 11_000_000),
        ];
        let maps = d.0.join("maps");
        assert_eq!(paint_chroms(&counts, &maps, 100, 5.0).unwrap(), ["chr20"]);
        let map = read_map(&maps.join("chr20.map")).unwrap();
        assert_eq!(cm_at(&map, 6_000_000), 5.0);
        assert_eq!(cm_at(&map, 1), 0.0);
        assert_eq!(cm_at(&map, 99_000_000), 10.0);
        assert_eq!(cm_at(&[], 5), 0.0);
        let c = d.put(
            "counts.tsv",
            "chr22\t500\t2000000\t8000000\nchr20\t500\t2000000\t8000000\nhomref\t3\n",
        );
        assert_eq!(
            run(&["chroms", "--counts", &c, "--maps", &maps.to_string_lossy()]),
            (0, "chr20 chr22\n".into(), String::new())
        );
        assert_eq!(
            run(&[
                "chroms",
                "--counts",
                &c,
                "--maps",
                &maps.to_string_lossy(),
                "--min-cm",
                "6.5"
            ])
            .1,
            "\n"
        );
    }
}
