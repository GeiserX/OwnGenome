//! `og family ibd`: the DNA segments two sequenced relatives share by
//! descent, for step 38 and the Nextflow FAMILY workflow
//! (`modules/local/ibd`).
//!
//! Reads, for each chromosome of `--chroms` (one per line): Beagle's phased
//! family VCF `PHASED/<chrom>.vcf.gz`, hap-ibd's segments
//! `IBD/<chrom>.ibd.gz` and the PLINK genetic map `MAPS/<chrom>.map`.
//! hap-ibd reports a segment per pair of haplotypes; here the segments of
//! one pair of people are merged into the stretches they share (one copy or
//! both). A phasing error breaks hap-ibd's segments, so two stretches of a
//! pair are joined when the markers between them hold at most
//! [`MAX_DISCORD`] opposite homozygote (one 0/0, the other 1/1), and a
//! pair's first and last stretch reach the chromosome's first and last
//! marker on the same test.
//!
//! Writes `ibd_segments.tsv` (one row per shared stretch),
//! `ibd_pairs.tsv` (one row per pair of members) and `ibd_report.md`.
//!
//! Each chromosome is reduced to its stretches as soon as it is read, so
//! memory holds one chromosome's markers at a time: a position and one
//! byte per member per marker. A missing or damaged input is an `ERROR:`
//! line and exit 2; each output is written whole or not at all.

use std::collections::HashMap;
use std::path::Path;

use og_io::{Lines, Utf8};
use og_pedigree::{Die, Kinship, Members, build_members, declared_relationship, read_pedigree};
use pycompat::fmt::fixed;
use pycompat::num::{parse_float, parse_int, sum as py_sum};
use pycompat::text::{split, strip};

use super::args::{Args, Kind, Spec};
use super::{Result, write_file};

pub const SPEC: Spec = Spec {
    prog: "og family ibd",
    usage: "usage: og family ibd [-h] --pedigree PEDIGREE --chroms CHROMS --maps MAPS --phased PHASED --ibd IBD --out-dir OUT_DIR\n",
    options: &[
        ("--pedigree", Kind::One),
        ("--chroms", Kind::One),
        ("--maps", Kind::One),
        ("--phased", Kind::One),
        ("--ibd", Kind::One),
        ("--out-dir", Kind::One),
    ],
    required: &[
        "--pedigree",
        "--chroms",
        "--maps",
        "--phased",
        "--ibd",
        "--out-dir",
    ],
};

const PAIR_COLUMNS: &str = "sample_a\tsample_b\tdeclared\texpected_relatedness\tsegments\tshared_cm\tibd2_cm\tlongest_cm\tanalysed_cm\tshared_fraction\trelatedness\tdegree\n";
const SEGMENT_COLUMNS: &str = "sample_a\tsample_b\tchrom\tstart\tend\tcm\tibd2_cm\n";
/// hap-ibd's default min-output: shorter segments are not reported.
const MIN_OUTPUT_CM: f64 = 2.0;
/// Opposite homozygotes allowed between two stretches that are joined.
pub const MAX_DISCORD: usize = 1;
/// Share of a first-degree pair's genome in both copies from which the
/// pair reads as full siblings (about 25%) rather than parent and child (0%).
pub const SIBLING_MIN_IBD2: f64 = 0.10;
const DEGREES: [(f64, &str); 5] = [
    (0.707, "same person or identical twins"),
    (0.354, "first degree"),
    (0.177, "second degree"),
    (0.0884, "third degree"),
    (0.0442, "fourth degree"),
];

/// A marker's ALT count: 0, 1 or 2; -1 when the genotype is not two 0/1 alleles.
type Dosage = i8;

/// `os.path.join(dir, name)` for a plain file name.
fn join(dir: &str, name: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') {
        format!("{dir}{name}")
    } else {
        format!("{dir}/{name}")
    }
}

/// The text lines of a file, plain or gzip, as Python reads them.
fn lines(path: &str) -> std::result::Result<Lines<Box<dyn std::io::BufRead>>, Die> {
    let f = og_io::open(path).map_err(|e| Die::io(path, &e))?;
    Ok(Lines::new(f, Utf8::Strict))
}

/// Calls `each(line number, line without its ending)` for every line.
fn for_lines(
    path: &str,
    mut each: impl FnMut(usize, &str) -> std::result::Result<(), Die>,
) -> std::result::Result<(), Die> {
    let mut l = lines(path)?;
    let mut line = String::new();
    let mut n = 0;
    while l.read_line(&mut line).map_err(|e| Die::io(path, &e))? {
        n += 1;
        each(n, line.trim_end_matches(['\r', '\n']))?;
    }
    Ok(())
}

/// Python's `bisect_left`, step for step.
fn bisect_left(a: &[i64], x: i64) -> usize {
    let (mut lo, mut hi) = (0, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < x { lo = mid + 1 } else { hi = mid }
    }
    lo
}

/// Python's `bisect_right`, step for step.
fn bisect_right(a: &[i64], x: i64) -> usize {
    let (mut lo, mut hi) = (0, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if x < a[mid] { hi = mid } else { lo = mid + 1 }
    }
    lo
}

/// A PLINK genetic map (chrom, id, cM, bp): cM at any position by linear
/// interpolation, held flat beyond the first and last row.
#[derive(Debug)]
pub struct GeneticMap {
    bp: Vec<i64>,
    cm: Vec<f64>,
}

impl GeneticMap {
    pub fn read(path: &str) -> std::result::Result<GeneticMap, Die> {
        let mut rows: Vec<(i64, f64)> = Vec::new();
        for_lines(path, |n, line| {
            let c = split(line, None, None);
            if c.is_empty() {
                return Ok(());
            }
            if c.len() < 4 {
                return Err(Die::one(format!(
                    "{path} line {n}: a PLINK map row has 4 columns (chrom, id, cM, bp)"
                )));
            }
            let bp = parse_int(c[3]).ok().and_then(|v| i64::try_from(v).ok());
            match (bp, parse_float(c[2])) {
                (Some(bp), Ok(cm)) => rows.push((bp, cm)),
                _ => {
                    return Err(Die::one(format!(
                        "{path} line {n}: cM and bp must be numbers"
                    )));
                }
            }
            Ok(())
        })?;
        GeneticMap::from_rows(path, rows)
    }

    fn from_rows(path: &str, mut rows: Vec<(i64, f64)>) -> std::result::Result<GeneticMap, Die> {
        if rows.len() < 2 {
            return Err(Die::one(format!(
                "{path}: a genetic map needs two rows or more"
            )));
        }
        // Python sorts the (bp, cM) tuples: by bp, then cM.
        rows.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        });
        let (bp, cm) = rows.into_iter().unzip();
        Ok(GeneticMap { bp, cm })
    }

    /// cM at `pos`.
    pub fn at(&self, pos: i64) -> f64 {
        let i = bisect_left(&self.bp, pos);
        if i == 0 {
            return self.cm[0];
        }
        if i == self.bp.len() {
            return self.cm[self.cm.len() - 1];
        }
        let (b0, b1) = (self.bp[i - 1], self.bp[i]);
        if b1 == b0 {
            return self.cm[i];
        }
        self.cm[i - 1] + (self.cm[i] - self.cm[i - 1]) * (pos - b0) as f64 / (b1 - b0) as f64
    }
}

/// The ALT count of a GT field: `0/0`, `0|1`, `1/1` and the like; -1 for
/// anything else (missing, multi-allelic, haploid).
pub fn dosage(field: &str) -> Dosage {
    let gt = field.split(':').next().unwrap_or("").as_bytes();
    match gt {
        [a @ (b'0' | b'1'), b'/' | b'|', b @ (b'0' | b'1')] => ((a - b'0') + (b - b'0')) as Dosage,
        _ => -1,
    }
}

/// One chromosome of the phased VCF: the samples, the marker positions and
/// each column's dosage per marker.
#[derive(Debug)]
pub struct Phased {
    pub samples: Vec<String>,
    pub pos: Vec<i64>,
    /// One vector per sample column, each as long as `pos`.
    pub dos: Vec<Vec<Dosage>>,
}

impl Phased {
    /// Reads `path`, the phased VCF of `chrom`; a record of another
    /// chromosome is an error.
    pub fn read(path: &str, chrom: &str) -> std::result::Result<Phased, Die> {
        let mut samples: Option<Vec<String>> = None;
        let mut pos = Vec::new();
        let mut dos: Vec<Vec<Dosage>> = Vec::new();
        for_lines(path, |n, line| {
            if line.starts_with("#CHROM") {
                if !pos.is_empty() {
                    return Err(Die::one(format!(
                        "{path} line {n}: a #CHROM line after the records"
                    )));
                }
                let s: Vec<String> = line.split('\t').skip(9).map(str::to_owned).collect();
                dos = vec![Vec::new(); s.len()];
                samples = Some(s);
            } else if !line.starts_with('#') {
                let Some(s) = &samples else {
                    return Err(Die::one(format!(
                        "{path} has no #CHROM line before its records"
                    )));
                };
                let mut c = line.split('\t');
                let here = c.next().unwrap_or_default();
                if here != chrom {
                    return Err(Die::one(format!(
                        "{path} line {n}: chromosome {here} in the file of {chrom}"
                    )));
                }
                let p = c
                    .next()
                    .and_then(|v| parse_int(v).ok())
                    .and_then(|v| i64::try_from(v).ok())
                    .ok_or_else(|| {
                        Die::one(format!("{path} line {n}: the position must be an integer"))
                    })?;
                pos.push(p);
                let mut k = 0;
                for (d, g) in dos.iter_mut().zip(c.skip(7)) {
                    d.push(dosage(g));
                    k += 1;
                }
                if k < s.len() {
                    return Err(Die::one(format!(
                        "{path} line {n}: {k} sample columns, the #CHROM line names {}",
                        s.len()
                    )));
                }
            }
            Ok(())
        })?;
        let Some(samples) = samples else {
            return Err(Die::one(format!("{path} has no #CHROM line")));
        };
        Ok(Phased { samples, pos, dos })
    }
}

/// One hap-ibd segment, its pair in name order: (a, b, hap of a, hap of b, start, end).
type IbdRow = (String, String, u8, u8, i64, i64);

/// hap-ibd's `.ibd.gz`: sample1, hap1, sample2, hap2, chrom, start, end, cM.
fn read_ibd(path: &str, chrom: &str) -> std::result::Result<Vec<IbdRow>, Die> {
    let mut rows = Vec::new();
    for_lines(path, |n, line| {
        if line.is_empty() {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 8 {
            return Err(Die::one(format!(
                "{path} line {n}: hap-ibd rows have 8 columns, this one {}",
                c.len()
            )));
        }
        if c[4] != chrom {
            return Err(Die::one(format!(
                "{path} line {n}: chromosome {} in the file of {chrom}",
                c[4]
            )));
        }
        let int = |v: &str| parse_int(v).ok().and_then(|v| i64::try_from(v).ok());
        let (Some(ha), Some(hb), Some(s), Some(e)) = (int(c[1]), int(c[3]), int(c[5]), int(c[6]))
        else {
            return Err(Die::one(format!(
                "{path} line {n}: haplotypes and positions must be integers"
            )));
        };
        if !(1..=2).contains(&ha) || !(1..=2).contains(&hb) || e < s {
            return Err(Die::one(format!(
                "{path} line {n}: a haplotype is 1 or 2 and the end is not before the start"
            )));
        }
        let (mut a, mut b, mut ha, mut hb) = (c[0], c[2], ha as u8, hb as u8);
        if a == b {
            return Ok(());
        }
        if b < a {
            (a, b, ha, hb) = (b, a, hb, ha);
        }
        rows.push((a.to_owned(), b.to_owned(), ha, hb, s, e));
        Ok(())
    })?;
    Ok(rows)
}

type Iv = (i64, i64);
/// One hap-ibd segment of a known pair: (hap of a, hap of b, start, end).
type HapSeg = (u8, u8, i64, i64);

/// The union of intervals: sorted, overlapping or touching ones merged.
pub fn union(mut iv: Vec<Iv>) -> Vec<Iv> {
    iv.sort_unstable();
    let mut out: Vec<Iv> = Vec::with_capacity(iv.len());
    for (s, e) in iv {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// The intersection of two sorted interval lists; empty overlaps dropped.
pub fn intersect(x: &[Iv], y: &[Iv]) -> Vec<Iv> {
    let (mut out, mut i, mut j) = (Vec::new(), 0, 0);
    while i < x.len() && j < y.len() {
        let (s, e) = (x[i].0.max(y[j].0), x[i].1.min(y[j].1));
        if s < e {
            out.push((s, e));
        }
        if x[i].1 < y[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}

/// Where both copies are shared: a's haplotype 1 with one of b's and a's
/// haplotype 2 with the other. `segs` holds (hap of a, hap of b, start, end).
pub fn ibd2(segs: &[HapSeg]) -> Vec<Iv> {
    let u = |ha: u8, hb: u8| {
        union(
            segs.iter()
                .filter(|r| (r.0, r.1) == (ha, hb))
                .map(|r| (r.2, r.3))
                .collect(),
        )
    };
    let mut both = intersect(&u(1, 1), &u(2, 2));
    both.extend(intersect(&u(1, 2), &u(2, 1)));
    union(both)
}

/// Markers strictly between positions `a` and `b` where one of the two is
/// 0/0 and the other 1/1.
fn opposite_homozygotes(pos: &[i64], dx: &[Dosage], dy: &[Dosage], a: i64, b: i64) -> usize {
    let (i, j) = (bisect_right(pos, a), bisect_left(pos, b));
    if i >= j {
        return 0;
    }
    dx[i..j]
        .iter()
        .zip(&dy[i..j])
        .filter(|&(&x, &y)| (x == 0 && y == 2) || (x == 2 && y == 0))
        .count()
}

/// Joins stretches with at most [`MAX_DISCORD`] opposite homozygotes
/// between them, and reaches the first and last marker on the same test.
pub fn stitch(iv: &[Iv], pos: &[i64], dx: &[Dosage], dy: &[Dosage]) -> Vec<Iv> {
    let Some((&first, rest)) = iv.split_first() else {
        return Vec::new();
    };
    let mut out = vec![first];
    for &(s, e) in rest {
        let last = out.last_mut().expect("one stretch");
        if opposite_homozygotes(pos, dx, dy, last.1, s) <= MAX_DISCORD {
            last.1 = last.1.max(e);
        } else {
            out.push((s, e));
        }
    }
    let (p0, p1) = (pos[0], pos[pos.len() - 1]);
    if opposite_homozygotes(pos, dx, dy, p0 - 1, out[0].0) <= MAX_DISCORD {
        out[0].0 = out[0].0.min(p0);
    }
    let n = out.len() - 1;
    if opposite_homozygotes(pos, dx, dy, out[n].1, p1 + 1) <= MAX_DISCORD {
        out[n].1 = out[n].1.max(p1);
    }
    out
}

/// The degree a relatedness points to.
fn degree(rel: f64, ibd2_frac: f64, n: usize) -> String {
    if n == 0 {
        return "no shared segment".to_owned();
    }
    for (lo, name) in DEGREES {
        if rel >= lo {
            if name == "first degree" {
                let kind = if ibd2_frac >= SIBLING_MIN_IBD2 {
                    " (full siblings)"
                } else {
                    " (parent-child)"
                };
                return format!("{name}{kind}");
            }
            return name.to_owned();
        }
    }
    "distant or unrelated (short segments only)".to_owned()
}

/// One shared stretch of a pair.
struct Stretch {
    chrom: usize,
    start: i64,
    end: i64,
    cm: f64,
    cm2: f64,
}

/// What the chromosomes give: their names, each one's analysed cM, the
/// members, and each pair's stretches in chromosome order.
struct Found {
    chroms: Vec<String>,
    glen: Vec<f64>,
    samples: Vec<String>,
    by_pair: HashMap<(String, String), Vec<Stretch>>,
}

/// Reads the chromosomes one at a time and keeps only their stretches.
fn find(a: &Args) -> std::result::Result<Found, Die> {
    let arg = |k: &str| a.one(k).unwrap_or_default().to_owned();
    let (chroms_path, maps, phased, ibd) = (
        arg("--chroms"),
        arg("--maps"),
        arg("--phased"),
        arg("--ibd"),
    );
    let mut chroms: Vec<String> = Vec::new();
    for_lines(&chroms_path, |n, line| {
        let c = strip(line, None);
        if c.is_empty() {
            return Ok(());
        }
        if chroms.iter().any(|x| x == c) {
            return Err(Die::one(format!(
                "{chroms_path} line {n}: {c} is listed twice"
            )));
        }
        chroms.push(c.to_owned());
        Ok(())
    })?;
    let mut found = Found {
        chroms: Vec::new(),
        glen: Vec::new(),
        samples: Vec::new(),
        by_pair: HashMap::new(),
    };
    let mut first: Option<Vec<String>> = None;
    for (ci, c) in chroms.iter().enumerate() {
        let m = GeneticMap::read(&join(&maps, &format!("{c}.map")))?;
        let ph = Phased::read(&join(&phased, &format!("{c}.vcf.gz")), c)?;
        match &first {
            None => first = Some(ph.samples.clone()),
            Some(s) if *s != ph.samples => {
                return Err(Die::one(format!(
                    "{c}.vcf.gz has the samples {}, {}.vcf.gz {}",
                    ph.samples.join(","),
                    chroms[0],
                    s.join(",")
                )));
            }
            Some(_) => {}
        }
        if ph.pos.is_empty() {
            return Err(Die::one(format!("{c}.vcf.gz has no marker")));
        }
        let (p0, p1) = (ph.pos[0], ph.pos[ph.pos.len() - 1]);
        found.glen.push(m.at(p1) - m.at(p0));
        // A name the VCF holds twice reads as its last column.
        let col: HashMap<&str, usize> = ph
            .samples
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        let mut segs: HashMap<(String, String), Vec<HapSeg>> = HashMap::new();
        let mut order: Vec<(String, String)> = Vec::new();
        for (x, y, ha, hb, s0, e0) in read_ibd(&join(&ibd, &format!("{c}.ibd.gz")), c)? {
            if !col.contains_key(x.as_str()) || !col.contains_key(y.as_str()) {
                return Err(Die::one(format!(
                    "{c}.ibd.gz names {x} and {y}; the phased VCF holds {}",
                    ph.samples.join(",")
                )));
            }
            let key = (x, y);
            if !segs.contains_key(&key) {
                order.push(key.clone());
            }
            segs.entry(key).or_default().push((ha, hb, s0, e0));
        }
        for key in order {
            let rows = &segs[&key];
            let (dx, dy) = (&ph.dos[col[key.0.as_str()]], &ph.dos[col[key.1.as_str()]]);
            let two = ibd2(rows);
            let all = union(rows.iter().map(|r| (r.2, r.3)).collect());
            let out = found.by_pair.entry(key).or_default();
            for (s0, e0) in stitch(&all, &ph.pos, dx, dy) {
                let cm = m.at(e0) - m.at(s0);
                let cm2 = py_sum(
                    intersect(&[(s0, e0)], &two)
                        .into_iter()
                        .map(|(s, e)| m.at(e) - m.at(s)),
                );
                out.push(Stretch {
                    chrom: ci,
                    start: s0,
                    end: e0,
                    cm,
                    cm2,
                });
            }
        }
    }
    found.samples = first.unwrap_or_default();
    found.samples.sort();
    found.chroms = chroms;
    Ok(found)
}

/// The relationship the pedigree declares for each pair, and the
/// relatedness it implies.
struct Declared {
    members: Members,
}

impl Declared {
    fn read(path: &str) -> std::result::Result<Declared, Die> {
        let (rows, mut errors) = read_pedigree(path)?;
        let (members, more) = build_members(&rows, &HashMap::new(), path, "the samplesheet");
        errors.extend(more);
        if !errors.is_empty() {
            return Err(Die::one(errors.join("; ")));
        }
        Ok(Declared { members })
    }

    fn pair(&self, kin: &mut Kinship<'_>, a: &str, b: &str) -> (String, String) {
        if !self.members.contains_key(a) || !self.members.contains_key(b) {
            return ("not in the pedigree".to_owned(), String::new());
        }
        let e = 2.0 * kin.phi(a, b);
        (declared_relationship(a, b, &self.members, e), fixed(e, 3))
    }
}

/// Runs `og family ibd`.
pub fn run(a: &Args) -> Result<()> {
    let found = find(a)?;
    let total = py_sum(found.glen.iter().copied());
    let declared = Declared::read(a.one("--pedigree").unwrap_or_default())?;
    let mut kin = Kinship::new(&declared.members);

    let mut seg_out = String::from(SEGMENT_COLUMNS);
    let mut pair_out = String::from(PAIR_COLUMNS);
    let mut table = String::new();
    let none: Vec<Stretch> = Vec::new();
    let s = &found.samples;
    for (i, x) in s.iter().enumerate() {
        for y in &s[i + 1..] {
            let (mut shared, mut both, mut longest, mut n) = (0.0f64, 0.0f64, 0.0f64, 0usize);
            for st in found.by_pair.get(&(x.clone(), y.clone())).unwrap_or(&none) {
                seg_out.push_str(&format!(
                    "{x}\t{y}\t{}\t{}\t{}\t{}\t{}\n",
                    found.chroms[st.chrom],
                    st.start,
                    st.end,
                    fixed(st.cm, 2),
                    fixed(st.cm2, 2)
                ));
                shared += st.cm;
                both += st.cm2;
                if st.cm > longest {
                    longest = st.cm;
                }
                n += 1;
            }
            let frac = if total > 0.0 { shared / total } else { 0.0 };
            let frac2 = if total > 0.0 { both / total } else { 0.0 };
            let rel = 2.0 * ((frac - frac2) / 4.0 + frac2 / 2.0);
            let (d, e) = declared.pair(&mut kin, x, y);
            let deg = degree(rel, frac2, n);
            pair_out.push_str(&format!(
                "{x}\t{y}\t{d}\t{e}\t{n}\t{}\t{}\t{}\t{}\t{}\t{}\t{deg}\n",
                fixed(shared, 2),
                fixed(both, 2),
                fixed(longest, 2),
                fixed(total, 2),
                fixed(frac, 4),
                fixed(rel, 3)
            ));
            table.push_str(&format!(
                "| {x} and {y} | {d} | {} cM | {}% | {} cM | {n} | {} cM | {deg} |\n",
                fixed(shared, 1),
                fixed(100.0 * frac, 1),
                fixed(longest, 1),
                fixed(both, 1)
            ));
        }
    }

    let mut md = String::from("# Shared DNA segments\n\n");
    let chroms = &found.chroms;
    if chroms.is_empty() {
        md.push_str("No chromosome was phased, so no shared segment was looked for.\n");
    } else {
        md.push_str(&format!(
            "hap-ibd looked for the stretches of DNA each pair of members shares by descent, on the joint \
             call phased by Beagle against the reference panel. Analysed: {} chromosome{} ({}), {} cM of \
             genetic map. Only stretches of {} cM or more count: shorter ones cannot be told apart from \
             chance sharing.\n",
            chroms.len(),
            if chroms.len() == 1 { "" } else { "s" },
            chroms.join(", "),
            fixed(total, 1),
            fixed(MIN_OUTPUT_CM, 0)
        ));
        if chroms.len() < 22 {
            md.push_str("\nFewer than the 22 autosomes were analysed, so the degree below is a rough estimate.\n");
        }
        md.push_str("\n| Pair | Declared | Shared | Of the analysed map | Longest | Stretches | Both copies | Estimate |\n|---|---|---|---|---|---|---|---|\n");
        md.push_str(&table);
        md.push_str(
            "\nEvery stretch is in ibd_segments.tsv (chromosome, start, end, cM, and the cM where both \
             copies are shared). Two stretches a phasing error split are joined when the markers between \
             them hold at most one opposite homozygote (one 0/0, the other 1/1), which two people who share \
             a copy never are bar a genotype error. A parent and child share one copy along the whole \
             genome; full siblings share both copies over about a quarter of it.\n",
        );
    }

    let out = Path::new(a.one("--out-dir").unwrap_or_default());
    std::fs::create_dir_all(out)
        .map_err(|e| Die::one(format!("cannot write {}: {e}", out.display())))?;
    write_file(&out.join("ibd_segments.tsv"), seg_out.as_bytes())?;
    write_file(&out.join("ibd_pairs.tsv"), pair_out.as_bytes())?;
    write_file(&out.join("ibd_report.md"), md.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests;
