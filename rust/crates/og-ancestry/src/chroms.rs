//! `chroms`: the chromosomes to phase and paint, those with N genotyped
//! SNVs or more (default 100) spanning X cM or more (default 5). A
//! chromosome with fewer genotyped SNVs or a shorter span is not painted.
//! The two defaults are a floor we chose, not a measured accuracy limit.
//! On a whole genome every autosome spans more than 40 cM. The map is read by
//! [`og_admix::map`], which the engines share.

use crate::{Args, Fail, Result, Spec};
use og_admix::map::GeneticMap;
use std::io::Write;
use std::path::Path;

pub(crate) const SPEC: Spec = Spec {
    verb: "chroms",
    usage: "usage: og ancestry-fine chroms --counts counts.tsv --maps DIR [--min-sites N] [--min-cm X]\n",
    options: &["--counts", "--maps", "--min-sites", "--min-cm"],
    required: &["--counts", "--maps"],
};

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
        let map = GeneticMap::read(&path).map_err(|e| Fail(e.0))?;
        if map.cm_at(*last) - map.cm_at(*first) >= min_cm {
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
        let map = GeneticMap::read(&maps.join("chr20.map")).unwrap();
        assert_eq!(map.cm_at(6_000_000), 5.0);
        assert_eq!(map.cm_at(1), 0.0);
        assert_eq!(map.cm_at(99_000_000), 10.0);
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
