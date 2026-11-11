//! `report`: the step's outputs, from ADMIXTURE's `.Q` and bootstrap
//! replicates of each level and RFMix's `.msp.tsv` of each chromosome.
//!
//! - `<S>_ancestry_fine.tsv`: each group's fraction with its bootstrap
//!   standard error and 2.5 to 97.5 percentile range, each continental
//!   group's share of the painted genome, a `left_out` row for each group
//!   the model leaves out and a `not_painted` row for each chromosome RFMix
//!   failed on (the paint folder's `failed.txt`).
//! - `<S>_painting.bed` and `<S>_painting.png`: the painting.
//! - `<S>_ancestry_fine.txt`: the plain text, also printed.
//!
//! Every input is read and checked first; the four files are then written
//! under temporary names and renamed into place.

use crate::panel::{left_out_groups, left_out_text};
use crate::{Args, COLOURS, ChromKey, Fail, Level, MIN_REFERENCES, Pending, Result, Spec, label};
use pycompat::fmt::{fixed, ljust, rjust};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) const SPEC: Spec = Spec {
    verb: "report",
    usage: "usage: og ancestry-fine report --sample S --model DIR --admixture DIR [--paint DIR] [--panel NAME] --out DIR\n",
    options: &[
        "--sample",
        "--model",
        "--admixture",
        "--paint",
        "--panel",
        "--out",
    ],
    required: &["--sample", "--model", "--admixture", "--out"],
};

/// GRCh38 lengths of chr1 to chr22, for the painting picture.
pub const CHROM_LEN: [u64; 22] = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
    138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
    83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
];

/// Grey where nothing is painted.
const UNPAINTED: [u8; 3] = [225, 225, 225];

/// One level of `levels.tsv`: its name, K and the group codes in the `.P`
/// column order.
#[derive(Debug, Clone)]
pub struct LevelRow {
    pub name: String,
    pub k: i64,
    pub codes: Vec<String>,
}

pub fn read_levels(model: &str) -> Result<Vec<LevelRow>> {
    let path = Path::new(model).join("levels.tsv");
    let mut out = Vec::new();
    let mut n = 0u64;
    crate::for_each_line(&path, |line| {
        n += 1;
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() >= 3 && !c[0].starts_with('#') {
            let k = crate::int(c[1]).ok_or_else(|| {
                Fail(format!(
                    "{}: line {n}: K is not a whole number: {}",
                    path.display(),
                    crate::repr(c[1])
                ))
            })?;
            out.push(LevelRow {
                name: c[0].to_owned(),
                k,
                codes: c[2].split(',').map(str::to_owned).collect(),
            });
        }
        Ok(())
    })?;
    if out.is_empty() {
        return Err(Fail(format!("{model}/levels.tsv lists no level")));
    }
    Ok(out)
}

/// The `q` percentile of `xs`, linear between the two nearest ranks.
pub fn percentile(xs: &[f64], q: f64) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    let mut xs = xs.to_vec();
    xs.sort_by(f64::total_cmp);
    let i = q * (xs.len() - 1) as f64;
    let lo = i.floor() as usize;
    let hi = (lo + 1).min(xs.len() - 1);
    Some(xs[lo] + (xs[hi] - xs[lo]) * (i - lo as f64))
}

/// One group's fraction at one level.
#[derive(Debug, Clone, PartialEq)]
pub struct Fraction {
    pub level: Level,
    pub code: String,
    pub fraction: f64,
    /// Bootstrap standard error and 2.5 and 97.5 percentiles, with two
    /// replicates or more.
    pub spread: Option<(f64, f64, f64)>,
}

/// The rows of one level, largest first, and the number of bootstrap
/// replicates.
pub fn fractions(
    level: Level,
    row: &LevelRow,
    adm_dir: &str,
    sample: &str,
) -> Result<(Vec<Fraction>, usize)> {
    let k = row.k;
    let name = &row.name;
    let q_name = format!("{sample}_{name}.{k}.Q");
    let q = crate::order::read_matrix(&Path::new(adm_dir).join(&q_name))?;
    if q.len() != 1 || q[0].len() as i64 != k {
        return Err(Fail(format!(
            "{adm_dir}/{q_name} is not one row of {k} values"
        )));
    }
    let boot_path = Path::new(adm_dir).join(format!("{sample}_{name}.boot.Q"));
    let boot = if boot_path.is_file() {
        crate::order::read_matrix(&boot_path)?
    } else {
        Vec::new()
    };
    if boot.iter().any(|r| r.len() as i64 != k) {
        return Err(Fail(format!(
            "{} has a row that is not {k} values",
            boot_path.display()
        )));
    }
    let mut rows: Vec<Fraction> = row
        .codes
        .iter()
        .zip(&q[0])
        .enumerate()
        .map(|(j, (code, &f))| {
            let col: Vec<f64> = boot.iter().map(|r| r[j]).collect();
            let spread = (col.len() >= 2).then(|| {
                let n = col.len() as f64;
                let mean = pycompat::num::sum(col.iter().copied()) / n;
                let var =
                    pycompat::num::sum(col.iter().map(|x| (x - mean) * (x - mean))) / (n - 1.0);
                (
                    var.sqrt(),
                    percentile(&col, 0.025).expect("two or more"),
                    percentile(&col, 0.975).expect("two or more"),
                )
            });
            Fraction {
                level,
                code: code.clone(),
                fraction: f,
                spread,
            }
        })
        .collect();
    rows.sort_by(|a, b| {
        b.fraction
            .partial_cmp(&a.fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok((rows, boot.len()))
}

/// One painted stretch of one haplotype: chromosome, 0-based start, end,
/// haplotype (1 or 2) and group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub hap: u8,
    pub code: String,
}

/// RFMix's `.msp.tsv`: the code of each group number, and the windows as
/// (chrom, start, end, haplotype 1's number, haplotype 2's number).
type Msp = (HashMap<i64, String>, Vec<(String, i64, i64, i64, i64)>);

fn read_msp(path: &Path) -> Result<Msp> {
    let mut codes = HashMap::new();
    let mut segs = Vec::new();
    let mut n = 0u64;
    let bad = |n: u64, v: &str| {
        Fail(format!(
            "{}: line {n}: not a whole number: {}",
            path.display(),
            crate::repr(v)
        ))
    };
    crate::for_each_line(path, |line| {
        n += 1;
        if let Some(rest) = line.strip_prefix("#Subpopulation order/codes:") {
            for item in crate::split_ws(rest) {
                let (name, num) = item.rsplit_once('=').unwrap_or(("", item));
                codes.insert(crate::int(num).ok_or_else(|| bad(n, num))?, name.to_owned());
            }
        } else if !line.starts_with('#') {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() >= 8 {
                let num = |v: &str| crate::int(v).ok_or_else(|| bad(n, v));
                segs.push((
                    c[0].to_owned(),
                    num(c[1])?,
                    num(c[2])?,
                    num(c[6])?,
                    num(c[7])?,
                ));
            }
        }
        Ok(())
    })?;
    if codes.is_empty() {
        return Err(Fail(format!(
            "{} has no '#Subpopulation order/codes' line",
            path.display()
        )));
    }
    Ok((codes, segs))
}

/// The painting: RFMix's windows, neighbouring windows of one haplotype
/// with the same group merged, sorted by chromosome, start and haplotype.
pub fn painting(paint_dir: &Path) -> Result<Vec<Segment>> {
    let mut names: Vec<String> = std::fs::read_dir(paint_dir)
        .map_err(|e| crate::cannot_read(paint_dir, &e))?
        .filter_map(|e| e.ok().and_then(|e| e.file_name().into_string().ok()))
        .filter(|n| n.ends_with(".msp.tsv"))
        .collect();
    names.sort_unstable();
    let mut rows: Vec<Segment> = Vec::new();
    for name in names {
        let path = paint_dir.join(&name);
        let (codes, segs) = read_msp(&path)?;
        for hap in [1u8, 2] {
            let mut cur: Option<Segment> = None;
            for (chrom, start, end, h1, h2) in &segs {
                let num = if hap == 1 { h1 } else { h2 };
                let code = codes.get(num).ok_or_else(|| {
                    Fail(format!(
                        "{}: group number {num} is not on the '#Subpopulation order/codes' line",
                        path.display()
                    ))
                })?;
                match &mut cur {
                    Some(c) if c.chrom == *chrom && c.code == *code && *start <= c.end + 1 => {
                        c.end = *end;
                    }
                    _ => {
                        rows.extend(cur.take());
                        cur = Some(Segment {
                            chrom: chrom.clone(),
                            start: start - 1,
                            end: *end,
                            hap,
                            code: code.clone(),
                        });
                    }
                }
            }
            rows.extend(cur);
        }
    }
    rows.sort_by(|a, b| {
        crate::chrom_key(&a.chrom)
            .cmp(&crate::chrom_key(&b.chrom))
            .then(a.start.cmp(&b.start))
            .then(a.hap.cmp(&b.hap))
    });
    Ok(rows)
}

/// The painted base pairs, and each group's share of them in order of
/// first appearance.
pub fn painting_shares(rows: &[Segment]) -> (i64, Vec<(String, f64)>) {
    let total: i64 = rows.iter().map(|r| r.end - r.start).sum();
    let mut by: Vec<(String, i64)> = Vec::new();
    for r in rows {
        match by.iter_mut().find(|(c, _)| *c == r.code) {
            Some((_, v)) => *v += r.end - r.start,
            None => by.push((r.code.clone(), r.end - r.start)),
        }
    }
    if total == 0 {
        return (0, Vec::new());
    }
    (
        total,
        by.into_iter()
            .map(|(c, v)| (c, v as f64 / total as f64))
            .collect(),
    )
}

/// A plain RGB PNG of `width` x `height` pixels, rows top to bottom.
pub fn png_bytes(width: u32, height: u32, rgb: &[u8]) -> Vec<u8> {
    use flate2::Compression;
    use flate2::write::ZlibEncoder;
    let stride = 3 * width as usize;
    let mut raw = Vec::with_capacity((stride + 1) * height as usize);
    for row in rgb.chunks_exact(stride) {
        raw.push(0); // filter: none
        raw.extend_from_slice(row);
    }
    let mut z = ZlibEncoder::new(Vec::new(), Compression::best());
    z.write_all(&raw).expect("writing to a Vec cannot fail");
    let idat = z.finish().expect("writing to a Vec cannot fail");
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.extend_from_slice(&[8, 2, 0, 0, 0]);
    let mut out = b"\x89PNG\r\n\x1a\n".to_vec();
    for (kind, data) in [
        (b"IHDR", &ihdr[..]),
        (b"IDAT", &idat[..]),
        (b"IEND", &[][..]),
    ] {
        out.extend_from_slice(&(data.len() as u32).to_be_bytes());
        let start = out.len();
        out.extend_from_slice(kind);
        out.extend_from_slice(data);
        let mut crc = flate2::Crc::new();
        crc.update(&out[start..]);
        out.extend_from_slice(&crc.sum().to_be_bytes());
    }
    out
}

/// Chromosomes 1 to 22 as two bars each (the two haplotypes), drawn to
/// scale, each painted segment in its group's colour (black for a group
/// without one), the rest grey.
pub fn painting_png(rows: &[Segment]) -> Vec<u8> {
    const WIDTH: usize = 1000;
    const BAR: usize = 6;
    const GAP: usize = 6;
    const STEP: usize = 2 * BAR + 1 + GAP;
    let height = GAP + 22 * STEP;
    let scale = (WIDTH - 20) as f64 / CHROM_LEN[0] as f64;
    let mut img = vec![255u8; WIDTH * height * 3];
    let mut fill = |y: usize, x0: usize, x1: usize, colour: [u8; 3]| {
        for yy in y..(y + BAR).min(height) {
            for x in x0..x1.min(WIDTH) {
                img[(yy * WIDTH + x) * 3..][..3].copy_from_slice(&colour);
            }
        }
    };
    for (i, &n) in CHROM_LEN.iter().enumerate() {
        for hap in 0..2 {
            let y = GAP + i * STEP + hap * (BAR + 1);
            fill(y, 10, 10 + (n as f64 * scale) as usize, UNPAINTED);
        }
    }
    for r in rows {
        let ChromKey::Number(c @ 1..=22) = crate::chrom_key(&r.chrom) else {
            continue;
        };
        let y = GAP + (c as usize - 1) * STEP + (usize::from(r.hap) - 1) * (BAR + 1);
        // Truncation toward zero, as Python's int().
        let x0 = 10 + (r.start as f64 * scale) as i64;
        let x1 = (x0 + 1).max(10 + (r.end as f64 * scale) as i64);
        let colour = COLOURS
            .iter()
            .find(|(code, _)| *code == r.code)
            .map_or([0, 0, 0], |&(_, c)| c);
        let clamp = |x: i64| x.clamp(0, (WIDTH - 10) as i64) as usize;
        fill(y, clamp(x0), clamp(x1), colour);
    }
    png_bytes(WIDTH as u32, height as u32, &img)
}

fn pct(x: Option<f64>) -> String {
    x.map_or_else(|| ".".to_owned(), |x| format!("{}%", fixed(100.0 * x, 0)))
}

fn fmt4(x: Option<f64>) -> String {
    x.map_or_else(|| ".".to_owned(), |x| fixed(x, 4))
}

const LIMITS: [&str; 6] = [
    "What this means: each name is a group of reference populations we put together and named,",
    "not a nationality. The fractions say which groups your genome is most similar to; the range",
    "after each is the 2.5 to 97.5 percentile of the bootstrap. Continental groups are well apart;",
    "neighbouring regions are not: splits within Iberia are weak for every method, and Greek and",
    "Southern Italian ancestry overlap (the panel has no Greek population at all). Read a regional",
    "fraction with its range, never as one exact percentage. See docs/26b-ancestry-fine.md.",
];

/// What the plain-text summary says.
pub struct Summary<'a> {
    pub sample: &'a str,
    pub panel: &'a str,
    pub fractions: &'a [Fraction],
    pub n_boot: &'a BTreeMap<Level, usize>,
    pub shares: &'a [(String, f64)],
    pub painted_bp: i64,
    pub n_sites: &'a str,
    pub painted_chroms: &'a [String],
    pub gone: &'a [(Level, String)],
    pub failed: &'a [String],
}

/// Shares largest first, ties in order of first appearance.
fn by_share(shares: &[(String, f64)]) -> Vec<&(String, f64)> {
    let mut v: Vec<&(String, f64)> = shares.iter().collect();
    v.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    v
}

pub fn text_block(s: &Summary<'_>) -> String {
    let mut out: Vec<String> = vec![
        format!("Fine-grained ancestry of {} (step 26b)", s.sample),
        format!("Reference panel: {}; {}", s.panel, s.n_sites),
        String::new(),
    ];
    if !s.gone.is_empty() {
        out.push(format!(
            "Left out of the model: {}. The panel has fewer than {MIN_REFERENCES} unrelated",
            left_out_text(s.gone)
        ));
        out.push("reference people for it, too few to stand for a population, so the model has no such group and".into());
        out.push("any ancestry from it is spread over the groups below.".into());
        out.push(String::new());
    }
    for level in Level::ALL {
        let rows: Vec<&Fraction> = s.fractions.iter().filter(|r| r.level == level).collect();
        if rows.is_empty() {
            continue;
        }
        let name = match level {
            Level::Continent => "Continental",
            Level::Region => "Regional",
        };
        out.push(format!(
            "{name} groups (ADMIXTURE projection on {} groups, {} bootstrap replicates):",
            rows.len(),
            s.n_boot.get(&level).copied().unwrap_or(0)
        ));
        for r in rows {
            let high = r.spread.map_or(0.0, |(_, _, h)| h);
            if level == Level::Region && r.fraction < 0.01 && high < 0.05 {
                continue;
            }
            let spread = r.spread.map_or(String::new(), |(_, lo, hi)| {
                format!("  ({} to {})", pct(Some(lo)), pct(Some(hi)))
            });
            out.push(format!(
                "  {} {}{spread}",
                ljust(label(level, &r.code), 50),
                rjust(&pct(Some(r.fraction)), 5)
            ));
        }
        out.push(String::new());
    }
    if !s.shares.is_empty() {
        out.push(format!(
            "Chromosome painting (RFMix, continental groups) of {}: share of the {} Mb painted on the two haplotypes",
            s.painted_chroms.join(", "),
            fixed(s.painted_bp as f64 / 1e6, 1)
        ));
        for (code, sh) in by_share(s.shares) {
            out.push(format!(
                "  {} {}",
                ljust(label(Level::Continent, code), 50),
                rjust(&pct(Some(*sh)), 5)
            ));
        }
    } else if s.failed.is_empty() {
        out.push(
            "Chromosome painting: none. No chromosome had enough genotyped SNVs spanning 5 cM."
                .into(),
        );
    } else {
        out.push("Chromosome painting: none.".into());
    }
    if !s.failed.is_empty() {
        out.push(format!(
            "  Not painted: RFMix failed on {} (the step's log has its last lines).",
            s.failed.join(", ")
        ));
    }
    out.push(String::new());
    out.extend(LIMITS.iter().map(|l| (*l).to_owned()));
    out.join("\n") + "\n"
}

/// Python's `os.path.normpath` on a POSIX path.
fn normpath(p: &str) -> String {
    if p.is_empty() {
        return ".".into();
    }
    let slashes = if p.starts_with("//") && !p.starts_with("///") {
        2
    } else {
        usize::from(p.starts_with('/'))
    };
    let mut comps: Vec<&str> = Vec::new();
    for c in p.split('/') {
        if c.is_empty() || c == "." {
            continue;
        }
        if c != ".." || (slashes == 0 && comps.is_empty()) || comps.last() == Some(&"..") {
            comps.push(c);
        } else if !comps.is_empty() {
            comps.pop();
        }
    }
    let out = format!("{}{}", "/".repeat(slashes), comps.join("/"));
    if out.is_empty() { ".".into() } else { out }
}

/// Python's `os.path.dirname` on a POSIX path.
fn dirname(p: &str) -> &str {
    let head = &p[..p.rfind('/').map_or(0, |i| i + 1)];
    if !head.is_empty() && head.bytes().any(|b| b != b'/') {
        head.trim_end_matches('/')
    } else {
        head
    }
}

/// The 0-based rows of a TSV after its header, split on tabs.
fn tsv_rows(path: &Path) -> Result<Vec<Vec<String>>> {
    Ok(crate::read_lines(path)?
        .iter()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect())
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let sample = a.req("--sample");
    let model = a.req("--model");
    let adm = a.req("--admixture");
    let paint = a.get("--paint").unwrap_or("");
    let panel = a.get("--panel").unwrap_or("");
    let levels = read_levels(model)?;
    let mut fr: Vec<Fraction> = Vec::new();
    let mut n_boot: BTreeMap<Level, usize> = BTreeMap::new();
    let mut codes: BTreeMap<Level, HashSet<String>> = BTreeMap::new();
    for row in &levels {
        let level = Level::parse(&row.name)
            .ok_or_else(|| Fail(format!("{model}/levels.tsv: unknown level {}", row.name)))?;
        let (rows, n) = fractions(level, row, adm, sample)?;
        fr.extend(rows);
        n_boot.insert(level, n);
        codes.insert(level, row.codes.iter().cloned().collect());
    }
    let rows = if !paint.is_empty() && Path::new(paint).is_dir() {
        painting(Path::new(paint))?
    } else {
        Vec::new()
    };
    let (painted_bp, shares) = painting_shares(&rows);
    let mut painted_chroms: Vec<String> = rows.iter().map(|r| r.chrom.clone()).collect();
    painted_chroms.sort_by(|a, b| crate::chrom_key(a).cmp(&crate::chrom_key(b)).then(a.cmp(b)));
    painted_chroms.dedup();
    let mut failed: Vec<String> = Vec::new();
    let failed_path = Path::new(paint).join("failed.txt");
    if !paint.is_empty() && failed_path.is_file() {
        failed = crate::read_lines(&failed_path)?
            .iter()
            .map(|l| pycompat::num::py_strip(l).to_owned())
            .filter(|l| !l.is_empty())
            .collect();
        failed.sort_by(|a, b| crate::chrom_key(a).cmp(&crate::chrom_key(b)).then(a.cmp(b)));
        failed.dedup();
    }
    let samples_path = Path::new(dirname(&normpath(model))).join("samples.tsv");
    let gone = if samples_path.is_file() {
        let samples: Vec<Vec<String>> = tsv_rows(&samples_path)?
            .into_iter()
            .filter(|r| r.len() >= 4)
            .collect();
        left_out_groups(&samples, &codes)
    } else {
        Vec::new()
    };
    let stats = Path::new(adm).join(format!("{sample}_sites.txt"));
    let n_sites = if stats.is_file() {
        let text = std::fs::read(&stats).map_err(|e| crate::cannot_read(&stats, &e))?;
        let text = String::from_utf8(text)
            .map_err(|_| Fail(format!("cannot read {}: not UTF-8", stats.display())))?;
        pycompat::num::py_strip(&text).to_owned()
    } else {
        ".".to_owned()
    };

    let mut table = String::from("level\tcode\tgroup\tfraction\tse\tlow\thigh\n");
    for r in &fr {
        let (se, lo, hi) = r
            .spread
            .map_or((None, None, None), |(s, l, h)| (Some(s), Some(l), Some(h)));
        table.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.level.name(),
            r.code,
            label(r.level, &r.code),
            fmt4(Some(r.fraction)),
            fmt4(se),
            fmt4(lo),
            fmt4(hi)
        ));
    }
    for (code, sh) in by_share(&shares) {
        table.push_str(&format!(
            "painting\t{code}\t{}\t{}\t.\t.\t.\n",
            label(Level::Continent, code),
            fmt4(Some(*sh))
        ));
    }
    for (level, code) in &gone {
        table.push_str(&format!(
            "left_out\t{code}\t{}\t.\t.\t.\t.\n",
            left_out_text(&[(*level, code.clone())])
        ));
    }
    for c in &failed {
        table.push_str(&format!("not_painted\t{c}\tRFMix failed\t.\t.\t.\t.\n"));
    }
    let mut bed = String::from("#chrom\tstart\tend\thaplotype\tgroup\n");
    for r in &rows {
        bed.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.chrom, r.start, r.end, r.hap, r.code
        ));
    }
    let text = text_block(&Summary {
        sample,
        panel,
        fractions: &fr,
        n_boot: &n_boot,
        shares: &shares,
        painted_bp,
        n_sites: &n_sites,
        painted_chroms: &painted_chroms,
        gone: &gone,
        failed: &failed,
    });
    let png = painting_png(&rows);

    let dir = PathBuf::from(a.req("--out"));
    std::fs::create_dir_all(&dir).map_err(|e| crate::cannot_write(&dir, &e))?;
    let base = |suffix: &str| dir.join(format!("{sample}{suffix}"));
    let mut pending = Pending::new();
    pending.put(base("_ancestry_fine.tsv"), table.as_bytes())?;
    pending.put(base("_painting.bed"), bed.as_bytes())?;
    pending.put(base("_painting.png"), &png)?;
    pending.put(base("_ancestry_fine.txt"), text.as_bytes())?;
    pending.commit_all()?;
    out.write_all(text.as_bytes())
        .map_err(|e| Fail(format!("cannot write stdout: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    #[test]
    fn percentiles_interpolate() {
        let xs = [3.0, 1.0, 2.0, 4.0];
        assert_eq!(percentile(&xs, 0.5), Some(2.5));
        assert_eq!(percentile(&xs, 0.0), Some(1.0));
        assert_eq!(percentile(&xs, 1.0), Some(4.0));
        assert_eq!(percentile(&[7.0], 0.975), Some(7.0));
        assert_eq!(percentile(&[], 0.5), None);
        // The 2.5 percentile of 0.80..0.90 by 0.01: rank 0.25.
        let col: Vec<f64> = (0..11).map(|i| 0.8 + f64::from(i) / 100.0).collect();
        assert_eq!(fixed(percentile(&col, 0.025).unwrap(), 4), "0.8025");
        assert_eq!(fixed(percentile(&col, 0.975).unwrap(), 4), "0.8975");
    }

    #[test]
    fn paths_as_python_reads_them() {
        assert_eq!(dirname(&normpath("model")), "");
        assert_eq!(dirname(&normpath("/a/b/model/")), "/a/b");
        assert_eq!(dirname(&normpath("p/./x/../model")), "p");
        assert_eq!(dirname(&normpath("/model")), "/");
        assert_eq!(normpath("../a//b/"), "../a/b");
        assert_eq!(normpath("//a"), "//a");
    }

    fn sample(d: &Dir) -> (String, String, String) {
        let model = d.path("model");
        d.put(
            "model/levels.tsv",
            "continent\t3\tEAS,WEA,AFR\nregion\t2\tIBE,ITA\n",
        );
        d.put(
            "samples.tsv",
            "sample\tpopulation\tcontinent\tregion\treference\tpaint\ni1\tIBS\tWEA\tIBE\t1\t1\nt1\tTSI\tWEA\tITA\t1\t1\np0\tPapuanHighlands\tOCE\tOCE\t0\t0\n",
        );
        d.put("admixture/S_continent.3.Q", "0.05 0.9 0.05\n");
        let boot: String = (0..11)
            .map(|i| {
                format!(
                    "0.05 {} {}\n",
                    fixed(0.8 + f64::from(i) / 100.0, 2),
                    fixed(0.15 - f64::from(i) / 100.0, 2)
                )
            })
            .collect();
        d.put("admixture/S_continent.boot.Q", &boot);
        d.put("admixture/S_region.2.Q", "0.6 0.4\n");
        d.put(
            "paint/chr20.msp.tsv",
            "#Subpopulation order/codes: AFR=0\tEAS=1\tWEA=2\n#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS.0\tS.1\n\
             chr20\t1000\t5000\t1\t2\t10\t2\t2\nchr20\t5000\t9000\t2\t3\t10\t2\t0\nchr20\t9000\t11000\t3\t4\t10\t2\t0\n",
        );
        d.put("paint/failed.txt", "chr19\n");
        (model, d.path("admixture"), d.path("paint"))
    }

    #[test]
    fn the_report() {
        let d = Dir::new("report");
        let (model, adm, paint) = sample(&d);
        let out = d.path("out");
        let (rc, stdout, err) = run(&[
            "report",
            "--sample",
            "S",
            "--model",
            &model,
            "--admixture",
            &adm,
            "--paint",
            &paint,
            "--panel",
            "test",
            "--out",
            &out,
        ]);
        assert_eq!((rc, err.as_str()), (0, ""));
        let table = d.read("out/S_ancestry_fine.tsv");
        let rows: Vec<Vec<&str>> = table.lines().map(|l| l.split('\t').collect()).collect();
        assert_eq!(
            rows[0],
            ["level", "code", "group", "fraction", "se", "low", "high"]
        );
        let cont: Vec<&Vec<&str>> = rows.iter().filter(|r| r[0] == "continent").collect();
        assert_eq!(
            cont.iter().map(|r| r[1]).collect::<Vec<_>>(),
            ["WEA", "EAS", "AFR"]
        );
        assert_eq!(cont[0][3..], ["0.9000", "0.0332", "0.8025", "0.8975"]);
        assert_eq!(
            rows.iter().find(|r| r[0] == "region").unwrap()[4..],
            [".", ".", "."]
        );
        let painted: Vec<(&str, &str)> = rows
            .iter()
            .filter(|r| r[0] == "painting")
            .map(|r| (r[1], r[3]))
            .collect();
        assert_eq!(painted, [("WEA", "0.7000"), ("AFR", "0.3000")]);
        let extra: Vec<&[&str]> = rows
            .iter()
            .filter(|r| r[0] == "left_out" || r[0] == "not_painted")
            .map(|r| &r[..3])
            .collect();
        assert_eq!(
            extra,
            [
                &["left_out", "OCE", "Oceania (continental group)"][..],
                &["left_out", "OCE", "Oceania (regional group)"][..],
                &["not_painted", "chr19", "RFMix failed"][..]
            ]
        );
        assert_eq!(
            d.read("out/S_painting.bed"),
            "#chrom\tstart\tend\thaplotype\tgroup\nchr20\t999\t11000\t1\tWEA\nchr20\t999\t5000\t2\tWEA\nchr20\t4999\t11000\t2\tAFR\n"
        );
        let text = d.read("out/S_ancestry_fine.txt");
        assert_eq!(stdout, text);
        assert!(
            text.contains("West Eurasia (Europe, Middle East, North Africa)     90%  (80% to 90%)")
        );
        assert!(text.contains("not a nationality") && text.contains("Greek and"));
        assert!(text.contains(
            "Left out of the model: Oceania (continental and regional group). The panel has fewer than 10 unrelated"
        ));
        assert!(text.contains("Not painted: RFMix failed on chr19") && text.contains("of chr20:"));
        let png = std::fs::read(d.0.join("out/S_painting.png")).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        assert_eq!(&png[16..24], &[0, 0, 3, 232, 0, 0, 1, 168]);
        // Without a painting the text says none.
        std::fs::remove_file(d.0.join("paint/failed.txt")).unwrap();
        std::fs::remove_file(d.0.join("paint/chr20.msp.tsv")).unwrap();
        let (rc, stdout, _) = run(&[
            "report",
            "--sample",
            "S",
            "--model",
            &model,
            "--admixture",
            &adm,
            "--paint",
            &paint,
            "--out",
            &out,
        ]);
        assert_eq!(rc, 0);
        assert!(stdout.contains("Chromosome painting: none. No chromosome had enough"));
    }

    #[test]
    fn the_painting_pixels() {
        let rows = vec![
            Segment {
                chrom: "chr1".into(),
                start: 0,
                end: 124_478_211,
                hap: 1,
                code: "AFR".into(),
            },
            Segment {
                chrom: "chr3".into(),
                start: 0,
                end: 1_000_000,
                hap: 1,
                code: "WEA".into(),
            },
            Segment {
                chrom: "chr2".into(),
                start: 0,
                end: 1,
                hap: 2,
                code: "ZZZ".into(),
            },
        ];
        let png = painting_png(&rows);
        let idat_at = png.windows(4).position(|w| w == b"IDAT").unwrap();
        let len = u32::from_be_bytes(png[idat_at - 4..idat_at].try_into().unwrap()) as usize;
        let mut raw = Vec::new();
        use std::io::Read;
        flate2::read::ZlibDecoder::new(&png[idat_at + 4..idat_at + 4 + len])
            .read_to_end(&mut raw)
            .unwrap();
        assert_eq!(raw.len(), 424 * 3001);
        let px = |x: usize, y: usize| &raw[y * 3001 + 1 + 3 * x..][..3];
        assert_eq!(px(0, 0), [255, 255, 255]);
        assert_eq!(px(10, 6), [230, 159, 0], "chr1 haplotype 1 starts African");
        assert_eq!(px(499, 6), [230, 159, 0]);
        assert_eq!(px(500, 6), [225, 225, 225], "and is grey past its half");
        assert_eq!(px(989, 6), [225, 225, 225]);
        assert_eq!(px(990, 6), [255, 255, 255]);
        assert_eq!(
            px(10, 6 + 19 + 7),
            [0, 0, 0],
            "a group without a colour is black, one pixel wide"
        );
        assert_eq!(px(11, 6 + 19 + 7), [225, 225, 225]);
        assert_eq!(
            px(10, 6 + 2 * 19),
            [0, 114, 178],
            "chr3 haplotype 1 starts West Eurasian"
        );
        assert_eq!(
            px(14, 6 + 2 * 19),
            [225, 225, 225],
            "for 1 Mb, under 4 pixels"
        );
    }

    #[test]
    fn a_bad_input_is_a_message() {
        let d = Dir::new("report-bad");
        let (model, adm, paint) = sample(&d);
        let out = d.path("out");
        let args = |adm: &str| {
            run(&[
                "report",
                "--sample",
                "S",
                "--model",
                &model,
                "--admixture",
                adm,
                "--paint",
                &paint,
                "--out",
                &out,
            ])
        };
        let (rc, _, err) = args(&d.path("none"));
        assert_eq!(rc, 1);
        assert!(
            err.starts_with("ERROR: cannot read ") && err.contains("S_continent.3.Q"),
            "{err}"
        );
        assert!(!d.0.join("out").exists(), "nothing is written");
        d.put("admixture/S_region.2.Q", "0.6 0.4 0.1\n");
        let (rc, _, err) = args(&adm);
        assert_eq!(rc, 1);
        assert!(
            err.ends_with("/S_region.2.Q is not one row of 2 values\n"),
            "{err}"
        );
        d.put("model/levels.tsv", "# nothing\n");
        assert!(args(&adm).2.ends_with("model/levels.tsv lists no level\n"));
        d.put("model/levels.tsv", "planet\t3\tA,B,C\n");
        assert!(
            args(&adm)
                .2
                .ends_with("model/levels.tsv: unknown level planet\n")
        );
    }
}
