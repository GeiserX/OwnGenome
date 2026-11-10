//! Runs of homozygosity: the `RG` rows of step 11's `bcftools roh` text.
//!
//! bcftools writes positions 1-based and inclusive; the map keeps spans
//! 0-based and end-exclusive, so `s = start - 1` and `e = end`. The `ST`
//! lines (one per site, most of a 200 MB file) are skipped by their first
//! two bytes, without decoding them.

use crate::facts::Sex;
use crate::reference::{Reference, chrom_index};
use crate::{Layer, State};
use std::io::BufRead;
use std::path::Path;

/// Where step 11 writes the text, under GENOME_DIR, in the order they are
/// tried: the single-step script, then the pipeline. Both are rows of
/// `assets/output_contract.tsv` (the contract-path test checks it).
pub const PATHS: [&str; 2] = [
    "{sample}/vcf/{sample}_roh.txt",
    "{sample}/roh/{sample}_roh.txt",
];

/// A run of 5 Mb or more with no artefact flag is a long run.
pub const LONG_BP: u64 = 5_000_000;

/// Why a run is not counted as long whatever its length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Art {
    None = 0,
    /// More than half of it lies in centromere or heterochromatin bands.
    Centromere = 1,
    /// On a male's single X, where every site reads homozygous.
    OneX = 2,
    /// On a male's single Y, for the same reason.
    OneY = 3,
}

/// `dict.art`, indexed by [`Art`].
pub const ART_NAMES: [&str; 4] = ["", "centromere", "one X", "one Y"];

#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    /// Index into [`crate::reference::CHROMS`].
    pub chrom: usize,
    /// 0-based start.
    pub s: u64,
    /// Exclusive end.
    pub e: u64,
    pub markers: u64,
    pub q: f64,
    pub long: bool,
    pub art: Art,
}

/// The layer and its rows, sorted by chromosome, start and end.
#[derive(Debug, Clone, PartialEq)]
pub struct Roh {
    pub layer: Layer,
    pub rows: Vec<Row>,
}

/// The input files this layer may read for `run`, under `root`, in order.
pub fn candidates(root: &Path, run: &str) -> Vec<(String, std::path::PathBuf)> {
    PATHS
        .iter()
        .map(|p| {
            let rel = p.replace("{sample}", run);
            let abs = root.join(&rel);
            (rel, abs)
        })
        .collect()
}

/// Reads the first of [`PATHS`] that exists for `run` under `root`.
pub fn read(root: &Path, run: &str, reference: &Reference, sex: Sex) -> Roh {
    let mut layer = Layer::new("11", "Runs of homozygosity");
    layer.caveat = crate::text::ROH_CAVEAT.to_owned();
    let all = candidates(root, run);
    let Some((rel, abs)) = all.iter().find(|(_, abs)| abs.is_file()) else {
        layer.source = all[0].0.clone();
        layer.state = State::Missing;
        layer.note = Some(crate::text::roh_missing(run));
        return Roh {
            layer,
            rows: Vec::new(),
        };
    };
    layer.source = rel.clone();
    let name = abs
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let parsed = og_io::open(abs)
        .map_err(|e| format!("{name}: {e}"))
        .and_then(|input| parse(input, &name, reference, sex));
    match parsed {
        Ok((rows, n_all)) => {
            layer.state = State::Ok;
            layer.count = rows.len();
            layer.totals = Some((n_all, rows.len()));
            layer.empty = Some(crate::text::ROH_EMPTY.to_owned());
            Roh { layer, rows }
        }
        Err(note) => {
            layer.state = State::Unreadable;
            layer.note = Some(note);
            Roh {
                layer,
                rows: Vec::new(),
            }
        }
    }
}

/// Parses a `bcftools roh` text. Returns the rows on the map's 24
/// chromosomes and the number of `RG` rows read, or `file:line: reason`.
pub fn parse(
    mut input: impl BufRead,
    name: &str,
    reference: &Reference,
    sex: Sex,
) -> Result<(Vec<Row>, usize), String> {
    let mut rows = Vec::new();
    let mut n_all = 0usize;
    let mut buf = Vec::with_capacity(256);
    let mut line_no = 0usize;
    loop {
        buf.clear();
        let n = input
            .read_until(b'\n', &mut buf)
            .map_err(|e| format!("{name}:{}: {e}", line_no + 1))?;
        if n == 0 {
            break;
        }
        line_no += 1;
        if !buf.starts_with(b"RG\t") {
            continue; // ST sites, comments and anything else bcftools adds
        }
        let at = |why: String| format!("{name}:{line_no}: {why}");
        let mut line = &buf[..];
        while let [rest @ .., b'\n' | b'\r'] = line {
            line = rest;
        }
        let line = std::str::from_utf8(line).map_err(|_| at("not UTF-8 text".to_owned()))?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 8 {
            return Err(at(format!("expected 8 columns, found {}", f.len())));
        }
        n_all += 1;
        let int = |i: usize, what: &str| {
            f[i].parse::<u64>()
                .map_err(|_| at(format!("{what} is not a whole number: {:?}", f[i])))
        };
        let (start, end, markers) = (int(3, "start")?, int(4, "end")?, int(6, "markers")?);
        if start == 0 || end < start {
            return Err(at(format!("{start}-{end} is not a 1-based span")));
        }
        let q = pycompat::num::parse_float(f[7])
            .ok()
            .filter(|q| q.is_finite())
            .ok_or_else(|| at(format!("quality is not a number: {:?}", f[7])))?;
        let Some(chrom) = chrom_index(f[2]) else {
            continue; // chrM, ALT and unplaced contigs are not on the map
        };
        let (s, e) = (start - 1, end);
        let c = &reference.chroms[chrom];
        let art = if sex == Sex::Male && c.name == "chrX" {
            Art::OneX
        } else if sex == Sex::Male && c.name == "chrY" {
            Art::OneY
        } else if 2 * c.masked_overlap(s, e) > e - s {
            Art::Centromere
        } else {
            Art::None
        };
        rows.push(Row {
            chrom,
            s,
            e,
            markers,
            q,
            long: e - s >= LONG_BP && art == Art::None,
            art,
        });
    }
    rows.sort_by_key(|r| (r.chrom, r.s, r.e));
    Ok((rows, n_all))
}

/// Bases covered per 1 Mb bin, per chromosome (index into `CHROMS`); only
/// chromosomes with a run are listed.
pub fn bins(rows: &[Row], reference: &Reference, size: u64) -> Vec<(usize, Vec<u64>)> {
    let mut out: Vec<(usize, Vec<u64>)> = Vec::new();
    for r in rows {
        if out.last().is_none_or(|(c, _)| *c != r.chrom) {
            let len = reference.chroms[r.chrom].len;
            out.push((r.chrom, vec![0; len.div_ceil(size) as usize]));
        }
        let bins = &mut out.last_mut().expect("pushed above").1;
        let mut a = r.s;
        while a < r.e {
            let k = (a / size) as usize;
            let b = r.e.min((k as u64 + 1) * size);
            if let Some(x) = bins.get_mut(k) {
                *x += b - a;
            }
            a = b;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference() -> Reference {
        Reference::load(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets")).unwrap()
    }

    const TEXT: &str = "\
# This file was produced by: bcftools roh
# ST\t[2]Sample\t[3]Chromosome\t[4]Position\t[5]State (0:HW, 1:AZ)\t[6]Quality
# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality
ST\tDEMO-001\tchr1\t10001\t0\t99.0
ST\tDEMO-001\tchr1\t10002\t1\t99.0\r
RG\tDEMO-001\tchr2\t1000001\t7000000\t6000000\t4100\t71.5
RG\tDEMO-001\tchr1\t1\t100\t100\t2\t30
RG\tDEMO-001\tchr1\t121000001\t123000000\t2000000\t900\t40.2\r
RG\tDEMO-001\tchrX\t20000001\t26000000\t6000000\t5000\t80
RG\tDEMO-001\tchrY\t3000001\t9000000\t6000000\t2100\t75
RG\tDEMO-001\tchrM\t1\t16569\t16569\t20\t50
RG\tDEMO-001\tchr13\t1\t20000000\t20000000\t9000\t60
";

    fn parse_text(t: &str, sex: Sex) -> Result<(Vec<Row>, usize), String> {
        parse(t.as_bytes(), "DEMO-001_roh.txt", &reference(), sex)
    }

    #[test]
    fn rg_rows_become_zero_based_spans() {
        let (rows, all) = parse_text(TEXT, Sex::Male).unwrap();
        assert_eq!(all, 7);
        let got: Vec<(&str, u64, u64, u64, bool, Art)> = rows
            .iter()
            .map(|r| {
                (
                    crate::reference::CHROMS[r.chrom],
                    r.s,
                    r.e,
                    r.markers,
                    r.long,
                    r.art,
                )
            })
            .collect();
        assert_eq!(
            got,
            [
                ("chr1", 0, 100, 2, false, Art::None),
                // 1.4 of its 2 Mb lie in chr1's centromere band.
                ("chr1", 121000000, 123000000, 900, false, Art::Centromere),
                ("chr2", 1000000, 7000000, 4100, true, Art::None),
                // 18.9 of its 20 Mb lie in chr13's short arm.
                ("chr13", 0, 20000000, 9000, false, Art::Centromere),
                ("chrX", 20000000, 26000000, 5000, false, Art::OneX),
                // A male's Y is hemizygous too, outside its masked spans.
                ("chrY", 3000000, 9000000, 2100, false, Art::OneY),
            ]
        );
        assert_eq!(rows[2].q, 71.5);
    }

    #[test]
    fn x_and_y_are_flagged_only_for_a_male() {
        for sex in [Sex::Female, Sex::Unknown] {
            let (rows, _) = parse_text(TEXT, sex).unwrap();
            for k in [22, 23] {
                let r = rows.iter().find(|r| r.chrom == k).unwrap();
                assert_eq!((r.art, r.long), (Art::None, true), "{sex:?} {k}");
            }
        }
    }

    #[test]
    fn bad_rows_name_file_and_line() {
        let t = "# c\nRG\tDEMO-001\tchr1\t5\t9\n";
        assert_eq!(
            parse_text(t, Sex::Unknown).unwrap_err(),
            "DEMO-001_roh.txt:2: expected 8 columns, found 5"
        );
        let t = "RG\tDEMO-001\tchr1\t9\t5\t1\t1\t1\n";
        assert_eq!(
            parse_text(t, Sex::Unknown).unwrap_err(),
            "DEMO-001_roh.txt:1: 9-5 is not a 1-based span"
        );
        // A site line that is not UTF-8 is skipped without being decoded.
        let t = b"ST\t\xff\nRG\tDEMO-001\tchr1\t1\t5\t5\tx\t1\n";
        assert_eq!(
            parse(&t[..], "r.txt", &reference(), Sex::Unknown).unwrap_err(),
            "r.txt:2: markers is not a whole number: \"x\""
        );
        let t = b"RG\tDEMO-001\tchr1\t1\t5\t5\t\xff\t1\n";
        assert_eq!(
            parse(&t[..], "r.txt", &reference(), Sex::Unknown).unwrap_err(),
            "r.txt:1: not UTF-8 text"
        );
    }

    #[test]
    fn bins_cover_each_base_once() {
        let r = reference();
        let rows = vec![Row {
            chrom: 0,
            s: 999_000,
            e: 2_000_500,
            markers: 1,
            q: 1.0,
            long: false,
            art: Art::None,
        }];
        let b = bins(&rows, &r, 1_000_000);
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].1.len(), 249);
        assert_eq!(&b[0].1[..4], [1000, 1_000_000, 500, 0]);
    }
}
