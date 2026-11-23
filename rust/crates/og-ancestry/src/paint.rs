//! `og ancestry-fine paint`: step 26b's chromosome painting, with our own
//! local-ancestry model ([`og_admix::paint`]) in place of RFMix.
//!
//! Inputs:
//! - `--phased DIR`: Beagle's phased chromosomes, every `<chrom>.vcf.gz`
//!   there, each holding the one person to paint;
//! - `--panel DIR` and `--maps DIR`: the panel's `<chrom>.vcf.gz` and the
//!   genetic map `<chrom>.map` (PLINK layout, the one Beagle reads) of each
//!   chromosome;
//! - `--model DIR`: `paint_map.tsv`, each painting reference and its
//!   continental group, and `levels.tsv`, whose continental level orders
//!   the groups;
//! - `--fractions DIR` with `--sample S`: the continental fractions of
//!   `og ancestry-fine project` (`<S>_continent.<K>.Q`), where the ancestry
//!   proportions of the model start. Without them every group starts
//!   equal.
//!
//! A chromosome is painted at the person's biallelic SNVs that the panel
//! holds with the same two alleles. Every chromosome is painted in one
//! model, whose parameters are fitted over the whole genome. A chromosome
//! that cannot be read (its phased file, its panel or its map missing or
//! broken, or no site shared with the panel) is named in `failed.txt`, with
//! the reason on stderr, and the others are still painted: `report` lists
//! it as not painted. If the model cannot be fitted, every chromosome goes
//! to `failed.txt`.
//!
//! Outputs in `--out-dir`, in RFMix's layout, so `report` reads them as it
//! read RFMix's:
//! - `<chrom>.msp.tsv`: one row for each run of sites where neither
//!   haplotype's group changes. The rows touch: the first starts at the
//!   first site, the last ends at the last site, and two neighbours meet
//!   halfway between the last site of one and the first site of the next,
//!   the second starting one base after the first ends. So `report` merges
//!   neighbouring rows of a haplotype with the same group into one stretch,
//!   with no gap between stretches.
//! - `params.tsv`: the fitted parameters
//!   ([`og_admix::paint::write_params`]).
//! - `failed.txt`: the chromosomes not painted, one per line, only when
//!   there is one.
//!
//! The files are written under temporary names and renamed into place once
//! the painting is done. The same inputs give the same bytes for any
//! `--threads`.

use crate::report::{LevelRow, read_levels, read_matrix};
use crate::{ArgError, Args, Fail, Pending, Result, Spec, repr};
use og_admix::map::GeneticMap;
use og_admix::paint::{
    Bits, ChromInput, DEFAULT_NE, DEFAULT_ROUNDS, Options, Params, RefPanel, write_params,
};
use og_admix::panel::Haplotypes;
use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) const SPEC: Spec = Spec {
    verb: "paint",
    usage: "usage: og ancestry-fine paint --phased DIR --panel DIR --maps DIR --model DIR --out-dir DIR [--fractions DIR --sample S] [--threads N]\n",
    options: &[
        "--phased",
        "--panel",
        "--maps",
        "--model",
        "--out-dir",
        "--fractions",
        "--sample",
        "--threads",
    ],
    required: &["--phased", "--panel", "--maps", "--model", "--out-dir"],
};

/// The least a group's start proportion is, as the fit keeps it, so no
/// group the projection put at 0 drops out of the model.
const MU_FLOOR: f64 = 1e-6;

/// `og ancestry-fine paint ARGS`: the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match Args::parse(&SPEC, args).and_then(|a| {
        a.int("--threads", 1)?;
        Ok(a)
    }) {
        Ok(a) => a,
        Err(ArgError::Help) => {
            let _ = write!(out, "{}", SPEC.usage);
            return 0;
        }
        Err(ArgError::Usage(msg)) => {
            let _ = write!(err, "{}", SPEC.usage);
            let _ = writeln!(err, "og ancestry-fine {}: error: {msg}", SPEC.verb);
            return 2;
        }
    };
    match main(&a, out, err) {
        Ok(()) => 0,
        Err(Fail(msg)) => {
            let _ = out.flush();
            let _ = writeln!(err, "ERROR: {msg}");
            1
        }
    }
}

fn fail(e: og_admix::Error) -> Fail {
    Fail(e.0)
}

/// `paint_map.tsv`: each painting reference and its group, the first two
/// columns of each line; blank lines and `#` lines skipped.
fn read_paint_map(path: &Path) -> Result<Vec<(String, String)>> {
    let mut rows = Vec::new();
    let mut seen = HashSet::new();
    let mut n = 0u64;
    crate::for_each_line(path, |line| {
        n += 1;
        if line.trim().is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = crate::split_ws(line).collect();
        if c.len() < 2 {
            return Err(Fail(format!(
                "{}: line {n}: a sample and its group, not {}",
                path.display(),
                repr(line)
            )));
        }
        if !seen.insert(c[0].to_owned()) {
            return Err(Fail(format!(
                "{}: line {n}: {} is named twice",
                path.display(),
                repr(c[0])
            )));
        }
        rows.push((c[0].to_owned(), c[1].to_owned()));
        Ok(())
    })?;
    if rows.is_empty() {
        return Err(Fail(format!("{}: no painting reference", path.display())));
    }
    Ok(rows)
}

/// The groups of the painting: those of the continental level that have a
/// painting reference, in the level's order, then any other group of the
/// references in order of first appearance.
fn groups(refs: &[(String, String)], continent: Option<&LevelRow>) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for code in continent.map_or(&[][..], |l| &l.codes[..]) {
        if refs.iter().any(|(_, g)| g == code) && !out.contains(code) {
            out.push(code.clone());
        }
    }
    for (_, g) in refs {
        if !out.contains(g) {
            out.push(g.clone());
        }
    }
    out
}

/// The start proportions: each group's continental fraction in
/// `<dir>/<sample>_continent.<K>.Q`, 0 for a group the level does not
/// have, kept at least [`MU_FLOOR`] and scaled to sum to 1.
fn start_mu(dir: &str, sample: &str, level: &LevelRow, groups: &[String]) -> Result<Vec<f64>> {
    let name = format!("{sample}_{}.{}.Q", level.name, level.k);
    let path = Path::new(dir).join(&name);
    let q = read_matrix(&path)?;
    if q.len() != 1 || q[0].len() as i64 != level.k || level.codes.len() as i64 != level.k {
        return Err(Fail(format!(
            "{dir}/{name} is not one row of {} values",
            level.k
        )));
    }
    if let Some(x) = q[0].iter().find(|x| !(x.is_finite() && **x >= 0.0)) {
        return Err(Fail(format!(
            "{dir}/{name} has a fraction of {x}; each must be 0 or more"
        )));
    }
    let mut mu: Vec<f64> = groups
        .iter()
        .map(|g| {
            level
                .codes
                .iter()
                .position(|c| c == g)
                .map_or(0.0, |j| q[0][j])
                .max(MU_FLOOR)
        })
        .collect();
    let total: f64 = mu.iter().sum();
    mu.iter_mut().for_each(|x| *x /= total);
    Ok(mu)
}

/// The sample names of a VCF's `#CHROM` line.
fn vcf_samples(path: &Path) -> Result<Vec<String>> {
    let input = og_io::open(path).map_err(|e| crate::cannot_read(path, &e))?;
    let mut lines = og_io::Lines::new(input, og_io::Utf8::Strict);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| crate::cannot_read(path, &e))?
    {
        let l = line.strip_suffix('\n').unwrap_or(&line);
        if l.starts_with("#CHROM") {
            return Ok(l.split('\t').skip(9).map(str::to_owned).collect());
        }
        if !l.starts_with('#') {
            break;
        }
    }
    Err(Fail(format!("{}: no #CHROM header line", path.display())))
}

/// One chromosome read for the painting, with each site's position and the
/// painted person's name.
struct Chrom {
    input: ChromInput,
    pos: Vec<u32>,
    sample: String,
}

/// Reads chromosome `c`: the person's two haplotypes at the sites the panel
/// holds with the same alleles, the references' haplotypes there, and each
/// site's cM.
fn read_chrom(
    c: &str,
    dirs: [&Path; 3],
    people: &[&str],
    person_group: &[usize],
    groups: usize,
) -> Result<Chrom> {
    let [phased, panel, maps] = dirs;
    let target_path = phased.join(format!("{c}.vcf.gz"));
    let samples = vcf_samples(&target_path)?;
    let [sample] = samples.as_slice() else {
        return Err(Fail(format!(
            "{} holds {} samples; the painting reads one person",
            target_path.display(),
            samples.len()
        )));
    };
    let target = Haplotypes::read(&target_path, &[sample.as_str()], None).map_err(fail)?;
    let sites = target.sites();
    if sites.is_empty() {
        return Err(Fail(format!(
            "{} has no biallelic SNV",
            target_path.display()
        )));
    }
    if target.chrom() != c {
        return Err(Fail(format!(
            "{} holds {}, not {c}",
            target_path.display(),
            target.chrom()
        )));
    }
    if let Some(w) = sites.windows(2).find(|w| w[1].pos <= w[0].pos) {
        return Err(Fail(format!(
            "{}: position {} after {}; the painting needs one SNV per position, in order",
            target_path.display(),
            w[1].pos,
            w[0].pos
        )));
    }
    let panel_path = panel.join(format!("{c}.vcf.gz"));
    let refs = Haplotypes::read(&panel_path, people, Some(sites)).map_err(fail)?;
    let shared = refs.sites();
    if shared.is_empty() {
        return Err(Fail(format!(
            "no SNV of {} is in {} with the same alleles",
            target_path.display(),
            panel_path.display()
        )));
    }
    let map = GeneticMap::read(&maps.join(format!("{c}.map"))).map_err(fail)?;
    let (mut h1, mut h2) = (Bits::default(), Bits::default());
    let mut cm = Vec::with_capacity(shared.len());
    let mut pos = Vec::with_capacity(shared.len());
    for s in shared {
        let m = sites
            .binary_search(s)
            .expect("the panel was read at the person's sites only");
        h1.push(target.allele(m, 0) == 1);
        h2.push(target.allele(m, 1) == 1);
        cm.push(map.cm_at(i64::from(s.pos)));
        pos.push(s.pos);
    }
    let panel = RefPanel::from_haplotypes(&refs, person_group, groups).map_err(fail)?;
    Ok(Chrom {
        input: ChromInput {
            name: c.to_owned(),
            cm,
            panel,
            targets: vec![h1, h2],
        },
        pos,
        sample: sample.clone(),
    })
}

/// The runs of sites where neither haplotype's call changes, as (first
/// site, last site).
fn runs(h1: &[usize], h2: &[usize]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut first = 0;
    for m in 1..=h1.len() {
        if m == h1.len() || h1[m] != h1[first] || h2[m] != h2[first] {
            out.push((first, m - 1));
            first = m;
        }
    }
    out
}

/// The base span of each run: the first starts at its first site, the last
/// ends at its last site, and two neighbours meet halfway between the last
/// site of one and the first site of the next, the second starting one
/// base after the first ends. Positions go up strictly, so each span holds
/// its sites.
fn spans(pos: &[u32], runs: &[(usize, usize)]) -> Vec<(u64, u64)> {
    let mut out = Vec::with_capacity(runs.len());
    let Some(&(first, _)) = runs.first() else {
        return out;
    };
    let mut start = u64::from(pos[first]);
    for (i, &(_, last)) in runs.iter().enumerate() {
        let end = match runs.get(i + 1) {
            Some(&(next, _)) => (u64::from(pos[last]) + u64::from(pos[next])) / 2,
            None => u64::from(pos[last]),
        };
        out.push((start, end));
        start = end + 1;
    }
    out
}

/// One chromosome's painting as an RFMix `.msp.tsv`: the group codes line,
/// the column header, then one row per run of sites where neither
/// haplotype's group changes, the rows touching (see [`spans`]). `sgpos`
/// and `egpos` are the cM of the row's first and last site, `n snps` its
/// sites, and the last two columns each haplotype's group number.
pub fn msp(
    chrom: &str,
    sample: &str,
    codes: &[String],
    pos: &[u32],
    cm: &[f64],
    calls: [&[usize]; 2],
) -> String {
    let order: Vec<String> = codes
        .iter()
        .enumerate()
        .map(|(i, c)| format!("{c}={i}"))
        .collect();
    let mut s = format!(
        "#Subpopulation order/codes: {}\n#chm\tspos\tepos\tsgpos\tegpos\tn snps\t{sample}.0\t{sample}.1\n",
        order.join("\t")
    );
    let rows = runs(calls[0], calls[1]);
    for (&(first, last), (spos, epos)) in rows.iter().zip(spans(pos, &rows)) {
        s.push_str(&format!(
            "{chrom}\t{spos}\t{epos}\t{:.2}\t{:.2}\t{}\t{}\t{}\n",
            cm[first],
            cm[last],
            last - first + 1,
            calls[0][first],
            calls[1][first]
        ));
    }
    s
}

/// The chromosomes of `dir`: every `<chrom>.vcf.gz`, in chromosome order.
fn phased_chroms(dir: &Path) -> Result<Vec<String>> {
    let entries = std::fs::read_dir(dir).map_err(|e| crate::cannot_read(dir, &e))?;
    let mut out = Vec::new();
    for e in entries {
        let e = e.map_err(|e| crate::cannot_read(dir, &e))?;
        if let Some(c) = e
            .file_name()
            .to_str()
            .and_then(|n| n.strip_suffix(".vcf.gz"))
            && !c.is_empty()
        {
            out.push(c.to_owned());
        }
    }
    out.sort_by(|a, b| crate::chrom_key(a).cmp(&crate::chrom_key(b)).then(a.cmp(b)));
    Ok(out)
}

fn main(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> Result<()> {
    let threads = a.int("--threads", 1).unwrap_or(1);
    if threads < 1 {
        return Err(Fail(format!(
            "--threads is {threads}; it must be 1 or more"
        )));
    }
    let fractions = match (a.get("--fractions"), a.get("--sample")) {
        (Some(d), Some(s)) => Some((d, s)),
        (None, None) => None,
        _ => return Err(Fail("--fractions and --sample go together".into())),
    };
    let model = a.req("--model");
    let refs = read_paint_map(&Path::new(model).join("paint_map.tsv"))?;
    let levels = read_levels(model)?;
    let continent = levels.iter().find(|l| l.name == "continent");
    let codes = groups(&refs, continent);
    let person_group: Vec<usize> = refs
        .iter()
        .map(|(_, g)| {
            codes
                .iter()
                .position(|c| c == g)
                .expect("every group listed")
        })
        .collect();
    let people: Vec<&str> = refs.iter().map(|(s, _)| s.as_str()).collect();
    let mu = match fractions {
        None => vec![1.0 / codes.len() as f64; codes.len()],
        Some((dir, sample)) => {
            let level = continent
                .ok_or_else(|| Fail(format!("{model}/levels.tsv has no continent level")))?;
            start_mu(dir, sample, level, &codes)?
        }
    };

    let phased = Path::new(a.req("--phased"));
    let dirs = [
        phased,
        Path::new(a.req("--panel")),
        Path::new(a.req("--maps")),
    ];
    let out_dir = PathBuf::from(a.req("--out-dir"));
    let chroms = phased_chroms(phased)?;
    std::fs::create_dir_all(&out_dir).map_err(|e| crate::cannot_write(&out_dir, &e))?;
    if chroms.is_empty() {
        let _ = writeln!(
            out,
            "Nothing to paint: no <chrom>.vcf.gz in {}",
            phased.display()
        );
        return Ok(());
    }

    let mut read: Vec<Chrom> = Vec::new();
    let mut failed: Vec<String> = Vec::new();
    for c in &chroms {
        match read_chrom(c, dirs, &people, &person_group, codes.len()) {
            Ok(ch) => read.push(ch),
            Err(Fail(why)) => {
                let _ = writeln!(err, "WARNING: {c} not painted: {why}");
                failed.push(c.clone());
            }
        }
    }

    let mut pending = Pending::new();
    if !read.is_empty() {
        let inputs: Vec<ChromInput> = read.iter().map(|c| c.input.clone()).collect();
        let sizes = inputs[0].panel.sizes().to_vec();
        let start = Params::initial(mu, &sizes, DEFAULT_NE);
        let opts = Options {
            threads: threads as usize,
            rounds: DEFAULT_ROUNDS,
        };
        match og_admix::paint::paint(&inputs, start, &opts) {
            Ok(painting) => {
                let groups: Vec<String> = codes
                    .iter()
                    .zip(&sizes)
                    .map(|(c, n)| format!("{c} {n}"))
                    .collect();
                let _ = writeln!(
                    out,
                    "Painted {} with {} reference haplotypes ({}): {} EM rounds, T {:.2} per Morgan, epsilon {:.4}",
                    read.iter()
                        .map(|c| format!("{} ({} SNVs)", c.input.name, c.pos.len()))
                        .collect::<Vec<_>>()
                        .join(", "),
                    sizes.iter().sum::<usize>(),
                    groups.join(", "),
                    painting.rounds,
                    painting.params.generations,
                    painting.params.epsilon
                );
                for (ch, p) in read.iter().zip(&painting.chroms) {
                    let text = msp(
                        &ch.input.name,
                        &ch.sample,
                        &codes,
                        &ch.pos,
                        &ch.input.cm,
                        [&p.haplotypes[0].calls, &p.haplotypes[1].calls],
                    );
                    pending.put(
                        out_dir.join(format!("{}.msp.tsv", ch.input.name)),
                        text.as_bytes(),
                    )?;
                }
                let mut params = Vec::new();
                write_params(&mut params, &painting, &codes).map_err(fail)?;
                pending.put(out_dir.join("params.tsv"), &params)?;
            }
            Err(e) => {
                let _ = writeln!(
                    err,
                    "WARNING: the painting model could not be fitted, so no chromosome is painted: {e}"
                );
                failed.extend(read.iter().map(|c| c.input.name.clone()));
            }
        }
    }
    if !failed.is_empty() {
        failed.sort_by(|a, b| crate::chrom_key(a).cmp(&crate::chrom_key(b)).then(a.cmp(b)));
        let text: String = failed.iter().map(|c| format!("{c}\n")).collect();
        pending.put(out_dir.join("failed.txt"), text.as_bytes())?;
    }
    pending.commit_all()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::painting;
    use crate::testutil::{Dir, run};
    use og_admix::rng::Rng;

    /// The sites of the fixed simulation, 1 kb apart from 1 kb on.
    const SITES: usize = 300;
    /// The site where the person's first haplotype turns from AFR to EAS.
    const SWITCH: usize = 150;

    /// Writes `text` gzipped at `name` under `d`.
    fn put_gz(d: &Dir, name: &str, text: &str) {
        let p = d.0.join(name);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        let f = std::fs::File::create(&p).unwrap();
        let mut e = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
    }

    fn vcf_header(samples: &[String]) -> String {
        format!(
            "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            samples.join("\t")
        )
    }

    /// A fixed two-group world on chr1: 10 AFR and 10 EAS painting
    /// references, one more panel person who paints nothing, 300 SNVs 1 kb
    /// apart at 0.1 cM per kb. At each site one group carries the ALT
    /// allele at 0.9 and the other at 0.1. The person T's first haplotype
    /// is drawn as AFR before site 150 and as EAS from it on, the second as
    /// EAS throughout. chr2 is phased but has no panel, so it cannot be
    /// painted.
    fn world(d: &Dir) -> Vec<String> {
        let mut rng = Rng::new(26);
        let mut people: Vec<String> = (0..10).map(|i| format!("a{i}")).collect();
        people.extend((0..10).map(|i| format!("e{i}")));
        people.push("x0".into());
        let mut panel = vcf_header(&people);
        let mut person = vcf_header(&["T".into()]);
        for m in 0..SITES {
            let p_afr = if rng.uniform() < 0.5 { 0.9 } else { 0.1 };
            let mut draw = |afr: bool| {
                let p = if afr { p_afr } else { 1.0 - p_afr };
                u8::from(rng.uniform() < p)
            };
            let pos = 1000 * (m + 1);
            panel.push_str(&format!("chr1\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT"));
            for i in 0..people.len() {
                let afr = i < 10;
                panel.push_str(&format!("\t{}|{}", draw(afr), draw(afr)));
            }
            panel.push('\n');
            person.push_str(&format!(
                "chr1\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT\t{}|{}\n",
                draw(m < SWITCH),
                draw(false)
            ));
        }
        put_gz(d, "panel/chr1.vcf.gz", &panel);
        put_gz(d, "phase/chr1.vcf.gz", &person);
        put_gz(
            d,
            "phase/chr2.vcf.gz",
            &(vcf_header(&["T".into()]) + "chr2\t5000\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\n"),
        );
        d.put(
            "maps/chr1.map",
            "chr1\t.\t0.0\t1000\nchr1\t.\t29.9\t300000\n",
        );
        d.put(
            "maps/chr2.map",
            "chr2\t.\t0.0\t1000\nchr2\t.\t10.0\t100000\n",
        );
        let map: String = people[..20]
            .iter()
            .enumerate()
            .map(|(i, p)| format!("{p}\t{}\n", if i < 10 { "AFR" } else { "EAS" }))
            .collect();
        d.put("model/paint_map.tsv", &map);
        d.put(
            "model/levels.tsv",
            "continent\t3\tEAS,WEA,AFR\nregion\t2\tWAF,EAS\n",
        );
        d.put("frc/T_continent.3.Q", "0.7 0.05 0.25\n");
        vec![
            "paint".into(),
            "--phased".into(),
            d.path("phase"),
            "--panel".into(),
            d.path("panel"),
            "--maps".into(),
            d.path("maps"),
            "--model".into(),
            d.path("model"),
            "--fractions".into(),
            d.path("frc"),
            "--sample".into(),
            "T".into(),
            "--out-dir".into(),
            d.path("paint"),
        ]
    }

    fn run_owned(args: &[String]) -> (u8, String, String) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        run(&args)
    }

    /// The painting of the fixed simulation, pinned: a change of the model,
    /// the fit or the writer shows here.
    const GOLDEN: &str = "#Subpopulation order/codes: EAS=0\tAFR=1\n\
                          #chm\tspos\tepos\tsgpos\tegpos\tn snps\tT.0\tT.1\n\
                          chr1\t1000\t150500\t0.00\t14.90\t150\t1\t0\n\
                          chr1\t150501\t300000\t15.00\t29.90\t150\t0\t0\n";

    #[test]
    fn a_fixed_simulation_paints_as_its_golden() {
        let d = Dir::new("paint-golden");
        let args = world(&d);
        let (rc, out, err) = run_owned(&args);
        assert_eq!(rc, 0, "{out}{err}");
        assert_eq!(d.read("paint/chr1.msp.tsv"), GOLDEN);
        // chr2 has no panel: named, and still not stopping chr1.
        assert_eq!(d.read("paint/failed.txt"), "chr2\n");
        assert!(
            err.contains("WARNING: chr2 not painted: cannot read "),
            "{err}"
        );
        assert!(
            out.contains("Painted chr1 (300 SNVs) with 40 reference haplotypes (EAS 20, AFR 20)"),
            "{out}"
        );
        let params = d.read("paint/params.tsv");
        assert!(
            params.starts_with("#parameter\tgroup\tvalue\nT\t.\t"),
            "{params}"
        );
        assert!(
            params.contains("\nmu\tEAS\t") && params.contains("\nem_rounds\t.\t"),
            "{params}"
        );
        // The same inputs on four threads give the same bytes.
        let d4 = Dir::new("paint-golden-4");
        let mut args = world(&d4);
        args.extend(["--threads".into(), "4".into()]);
        assert_eq!(run_owned(&args).0, 0);
        assert_eq!(d4.read("paint/chr1.msp.tsv"), d.read("paint/chr1.msp.tsv"));
        assert_eq!(d4.read("paint/params.tsv"), params);
    }

    #[test]
    fn the_report_finds_the_planted_breakpoint() {
        let d = Dir::new("paint-truth");
        let args = world(&d);
        assert_eq!(run_owned(&args).0, 0);
        let rows = painting(&d.0.join("paint")).unwrap();
        let hap = |h: u8| -> Vec<(i64, i64, &str)> {
            rows.iter()
                .filter(|r| r.hap == h)
                .map(|r| (r.start, r.end, r.code.as_str()))
                .collect()
        };
        // The second haplotype is EAS from the first site to the last.
        assert_eq!(hap(2), [(999, 300_000, "EAS")]);
        // The first turns from AFR to EAS within 5 kb (0.5 cM) of the
        // truth, halfway between sites 150 and 151 (150,500).
        let h1 = hap(1);
        assert_eq!(h1.len(), 2, "{h1:?}");
        assert_eq!(
            (h1[0].0, h1[0].2, h1[1].1, h1[1].2),
            (999, "AFR", 300_000, "EAS")
        );
        assert_eq!(h1[0].1, h1[1].0);
        assert!((h1[0].1 - 150_500).abs() <= 5_000, "{h1:?}");
    }

    /// Whether the report's painting of `msp` tiles each haplotype from
    /// 999 to 400 with no gap or overlap, and merges each haplotype's run
    /// of one group into one stretch.
    fn tiles(msp: &str) -> bool {
        let d = Dir::new("paint-tiles");
        d.put("chr1.msp.tsv", msp);
        let rows = painting(&d.0).unwrap();
        [1u8, 2].iter().all(|&h| {
            let segs: Vec<_> = rows.iter().filter(|r| r.hap == h).collect();
            let touching = segs
                .windows(2)
                .all(|w| w[0].end == w[1].start && w[0].code != w[1].code);
            touching
                && segs.first().is_some_and(|s| s.start == 99)
                && segs.last().is_some_and(|s| s.end == 400)
        })
    }

    /// Four sites; the first haplotype turns from group 0 to group 1 at the
    /// third, the second stays group 1: two rows.
    fn four_sites() -> (Vec<u32>, Vec<f64>, Vec<usize>, Vec<usize>) {
        (
            vec![100, 200, 300, 400],
            vec![0.1, 0.2, 0.3, 0.4],
            vec![0, 0, 1, 1],
            vec![1, 1, 1, 1],
        )
    }

    #[test]
    fn rows_touch_so_the_report_merges_them() {
        let (pos, cm, h1, h2) = four_sites();
        let codes = ["AFR".to_owned(), "EAS".to_owned()];
        let ours = msp("chr1", "T", &codes, &pos, &cm, [&h1, &h2]);
        assert_eq!(
            ours,
            "#Subpopulation order/codes: AFR=0\tEAS=1\n#chm\tspos\tepos\tsgpos\tegpos\tn snps\tT.0\tT.1\n\
             chr1\t100\t250\t0.10\t0.20\t2\t0\t1\nchr1\t251\t400\t0.30\t0.40\t2\t1\t1\n"
        );
        assert!(tiles(&ours));
        // Red control: a writer that ends each row at its last site and
        // starts the next at its first leaves a gap of 99 bases, and the
        // report then keeps the second haplotype's two rows apart.
        let gapped = ours
            .replace("chr1\t100\t250", "chr1\t100\t200")
            .replace("chr1\t251\t400", "chr1\t300\t400");
        assert_ne!(gapped, ours);
        assert!(!tiles(&gapped));
    }

    #[test]
    fn one_site_and_a_change_at_every_site() {
        assert_eq!(runs(&[3], &[1]), [(0, 0)]);
        assert_eq!(spans(&[7], &[(0, 0)]), [(7, 7)]);
        let r = runs(&[0, 1, 0], &[0, 0, 0]);
        assert_eq!(r, [(0, 0), (1, 1), (2, 2)]);
        // Neighbouring positions meet at the first: 10-10, 11-11, then 12-20.
        assert_eq!(spans(&[10, 11, 20], &r), [(10, 10), (11, 15), (16, 20)]);
    }

    #[test]
    fn no_phased_chromosome_paints_nothing() {
        let d = Dir::new("paint-none");
        let mut args = world(&d);
        std::fs::remove_dir_all(d.0.join("phase")).unwrap();
        std::fs::create_dir_all(d.0.join("phase")).unwrap();
        let (rc, out, _) = run_owned(&args);
        assert_eq!(rc, 0);
        assert!(out.starts_with("Nothing to paint"), "{out}");
        assert_eq!(std::fs::read_dir(d.0.join("paint")).unwrap().count(), 0);
        // --fractions needs --sample.
        let i = args.iter().position(|a| a == "--sample").unwrap();
        args.drain(i..i + 2);
        let (rc, _, err) = run_owned(&args);
        assert_eq!(rc, 1);
        assert!(
            err.contains("ERROR: --fractions and --sample go together"),
            "{err}"
        );
    }

    #[test]
    fn the_start_takes_the_projection_and_keeps_every_group() {
        let level = LevelRow {
            name: "continent".into(),
            k: 3,
            codes: vec!["EAS".into(), "WEA".into(), "AFR".into()],
        };
        let d = Dir::new("paint-mu");
        d.put("S_continent.3.Q", "0.6 0 0.4\n");
        let mu = start_mu(
            &d.0.to_string_lossy(),
            "S",
            &level,
            &["AFR".into(), "WEA".into(), "OCE".into()],
        )
        .unwrap();
        let total = 0.4 + MU_FLOOR + MU_FLOOR;
        assert_eq!(mu, [0.4 / total, MU_FLOOR / total, MU_FLOOR / total]);
        d.put("S_continent.3.Q", "0.6 0.4\n");
        let e = start_mu(&d.0.to_string_lossy(), "S", &level, &["AFR".into()]).unwrap_err();
        assert!(
            e.0.ends_with("S_continent.3.Q is not one row of 3 values"),
            "{}",
            e.0
        );
    }
}
