//! `og ancestry-fine project`: a person's ancestry fractions projected on
//! the groups' frequencies of the model, with a block bootstrap of the
//! genome and of the reference people. It replaces `admixture -P`; the
//! model and the solver are described in [`og_admix::project`].
//!
//! Inputs: the genotypes step 26b writes from the gVCF (`--vcf`), and the
//! model folder setup.sh writes (`--model`): `ref.bim`, `levels.tsv` and
//! each level's `<level>.P` (the frequency of the `.bim`'s A2, the VCF REF,
//! per group); for the bootstrap also `ref.geno`, `ref.fam`, the panel's
//! `samples.tsv` beside the model folder, and the genetic maps (`--maps`).
//!
//! A SNV counts when the person has a called genotype of the `.bim`'s two
//! alleles there; the last such record at a position decides it. With
//! fewer than `--min-sites` (100) the command stops.
//!
//! The bootstrap cuts the genotyped SNVs into blocks of 5 cM on the genetic
//! map, per chromosome, so no block holds two chromosomes. Each replicate
//! draws as many blocks as there are, with replacement, draws each group's
//! reference people again from `ref.geno`, counts the frequencies again and
//! projects again. Replicate r draws from its own stream of `--seed`, so
//! the same seed gives the same bytes whatever `--threads` is. A model
//! without `ref.geno` (made before og counted the frequencies) has no
//! people to draw: the bootstrap stops and says to run setup again.
//!
//! In `--out-dir`, per level: `<S>_<level>.<K>.Q` (one row of K fractions
//! in `levels.tsv` order, ADMIXTURE's layout), `<S>_<level>.boot.Q` (one row
//! per replicate) and `<S>_<level>.fit.tsv` (how the fit went); and
//! `<S>_sites.txt`, the SNVs genotyped. Everything is written under
//! temporary names and renamed into place once every level and replicate
//! is done.

use crate::freqs::a2_freq;
use crate::report::{LevelRow, read_levels};
use crate::{ArgError, Args, Fail, Level, Pending, Result, Spec, repr};
use og_admix::geno::RefGeno;
use og_admix::map::GeneticMap;
use og_admix::project::{self, BLOCK_CM, Options, Problem};
use og_admix::rng::Rng;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) const SPEC: Spec = Spec {
    verb: "project",
    usage: "usage: og ancestry-fine project --model DIR --vcf genotypes.vcf.gz --sample S --out-dir DIR [--maps DIR] [--bootstrap B] [--seed N] [--threads N] [--min-sites N] [--max-iterations N]\n",
    options: &[
        "--model",
        "--vcf",
        "--sample",
        "--out-dir",
        "--maps",
        "--bootstrap",
        "--seed",
        "--threads",
        "--min-sites",
        "--max-iterations",
    ],
    required: &["--model", "--vcf", "--sample", "--out-dir"],
};

/// `og ancestry-fine project ARGS`: the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match Args::parse(&SPEC, args).and_then(|a| {
        for name in [
            "--bootstrap",
            "--seed",
            "--threads",
            "--min-sites",
            "--max-iterations",
        ] {
            a.int(name, 0)?;
        }
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
    match main(&a, out) {
        Ok(()) => 0,
        Err(Fail(msg)) => {
            let _ = out.flush();
            let _ = writeln!(err, "ERROR: {msg}");
            1
        }
    }
}

/// The reference `.bim`: its text in one buffer, and each SNV as offsets
/// into it, so a panel of any size costs a few dozen bytes per SNV.
struct Bim {
    text: String,
    rows: Vec<BimRow>,
}

/// One SNV: the chromosome without `chr`, the alleles A1 and A2 as spans
/// of [`Bim::text`], and the position.
#[derive(Debug, Clone, Copy)]
struct BimRow {
    chrom: (u32, u32),
    a1: (u32, u32),
    a2: (u32, u32),
    pos: i64,
}

impl Bim {
    fn at(&self, (a, b): (u32, u32)) -> &str {
        &self.text[a as usize..b as usize]
    }
}

/// plink2 writes chr22 as 22 in a `.bim`; the VCFs say chr22.
fn norm_chrom(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

fn read_bim(path: &Path) -> Result<Bim> {
    let mut bim = Bim {
        text: String::new(),
        rows: Vec::new(),
    };
    let mut n = 0u64;
    let too_big = || Fail(format!("{}: larger than 4 GB", path.display()));
    crate::for_each_line(path, |line| {
        n += 1;
        let c: Vec<&str> = crate::split_ws(line).collect();
        if c.len() < 6 {
            return Ok(());
        }
        let pos = crate::int(c[3]).ok_or_else(|| {
            Fail(format!(
                "{}: line {n}: position is not a whole number: {}",
                path.display(),
                repr(c[3])
            ))
        })?;
        let mut span = |s: &str| -> Result<(u32, u32)> {
            let a = bim.text.len();
            bim.text.push_str(s);
            Ok((
                u32::try_from(a).map_err(|_| too_big())?,
                u32::try_from(bim.text.len()).map_err(|_| too_big())?,
            ))
        };
        let row = BimRow {
            chrom: span(norm_chrom(c[0]))?,
            a1: span(c[4])?,
            a2: span(c[5])?,
            pos,
        };
        bim.rows.push(row);
        Ok(())
    })?;
    Ok(bim)
}

/// A GT of exactly two called alleles, phased or not, as written.
fn gt_alleles(gt: &str) -> Option<[&str; 2]> {
    let mut parts = gt.split(['/', '|']);
    let (a, b) = (parts.next()?, parts.next()?);
    if parts.next().is_some() || a == "." || b == "." {
        return None;
    }
    Some([a, b])
}

/// The copies of A2 in a genotype (the record's alleles and the two
/// allele indices as written) for a `.bim` row with alleles A1, A2; `None`
/// for no call: an index that is no allele, or an allele the `.bim` does
/// not have.
pub fn a2_copies(alleles: &[&str], idx: [&str; 2], a1: &str, a2: &str) -> Option<u8> {
    let (mut n1, mut n2) = (0u8, 0u8);
    for i in idx {
        let i = crate::int(i)?;
        let called = *alleles.get(usize::try_from(i).ok()?)?;
        n1 += u8::from(called == a1);
        n2 += u8::from(called == a2);
    }
    (n1 + n2 == 2).then_some(n2)
}

/// The A2 copies at each `.bim` row, streamed from the VCF: the last
/// record with a called GT at a row's position decides it.
fn genotypes_at(bim: &Bim, vcf: &Path) -> Result<Vec<Option<u8>>> {
    // Every row by (chromosome, position), sorted, so a record finds its
    // rows by binary search (a .bim may list two SNVs at one position).
    let mut names: HashMap<&str, u32> = HashMap::new();
    let mut index: Vec<(u32, i64, u32)> = Vec::with_capacity(bim.rows.len());
    for (i, r) in bim.rows.iter().enumerate() {
        let next = u32::try_from(names.len()).expect("fewer chromosomes than rows");
        let c = *names.entry(bim.at(r.chrom)).or_insert(next);
        let i = u32::try_from(i).map_err(|_| Fail("more than 4 billion SNVs".into()))?;
        index.push((c, r.pos, i));
    }
    index.sort_unstable();
    let mut out: Vec<Option<u8>> = vec![None; bim.rows.len()];
    let mut n = 0u64;
    crate::for_each_line(vcf, |line| {
        n += 1;
        if line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 10 {
            return Ok(());
        }
        let Some(gt_at) = c[8].split(':').position(|f| f == "GT") else {
            return Ok(());
        };
        let Some(idx) = c[9].split(':').nth(gt_at).and_then(gt_alleles) else {
            return Ok(());
        };
        let pos = crate::int(c[1]).ok_or_else(|| {
            Fail(format!(
                "{}: line {n}: POS is not a whole number: {}",
                vcf.display(),
                repr(c[1])
            ))
        })?;
        let Some(&chrom) = names.get(norm_chrom(c[0])) else {
            return Ok(());
        };
        let from = index.partition_point(|&(ch, p, _)| (ch, p) < (chrom, pos));
        let rows = index[from..]
            .iter()
            .take_while(|&&(ch, p, _)| (ch, p) == (chrom, pos));
        let mut alleles = vec![c[3]];
        if c[4] != "." {
            alleles.extend(c[4].split(','));
        }
        for &(_, _, i) in rows {
            let r = &bim.rows[i as usize];
            out[i as usize] = a2_copies(&alleles, idx, bim.at(r.a1), bim.at(r.a2));
        }
        Ok(())
    })?;
    Ok(out)
}

/// One level's frequencies at the kept SNVs (row-major, SNV by group),
/// streamed from `<level>.P`; it must have one row of K numbers per `.bim`
/// SNV.
fn read_p(path: &Path, k: usize, keep: &[usize], bim_rows: usize) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(keep.len() * k);
    let (mut rows, mut next) = (0usize, keep.iter().peekable());
    crate::for_each_line(path, |line| {
        if next.peek() == Some(&&rows) {
            next.next();
            let before = out.len();
            for x in crate::split_ws(line) {
                out.push(pycompat::num::parse_float(x).map_err(|_| {
                    Fail(format!(
                        "{}: line {}: {} is not a number",
                        path.display(),
                        rows + 1,
                        repr(x)
                    ))
                })?);
            }
            if out.len() - before != k {
                return Err(Fail(format!(
                    "{}: line {}: {} values, not {k}",
                    path.display(),
                    rows + 1,
                    out.len() - before
                )));
            }
        }
        rows += 1;
        Ok(())
    })?;
    if rows != bim_rows {
        return Err(Fail(format!(
            "{} has {rows} rows but ref.bim {bim_rows} SNVs",
            path.display()
        )));
    }
    Ok(out)
}

/// A level of the model: its row of `levels.tsv` and its frequencies.
struct LevelModel {
    row: LevelRow,
    level: Level,
    k: usize,
    freqs: Vec<f64>,
    /// The people of each group, as `ref.geno` numbers them (bootstrap only).
    groups: Vec<Vec<usize>>,
}

/// The group of each `ref.fam` person at each level, from the panel's
/// `samples.tsv` beside the model folder.
fn assign_groups(model: &Path, levels: &mut [LevelModel], people: usize) -> Result<()> {
    let fam_path = model.join("ref.fam");
    let mut ids = Vec::new();
    crate::for_each_line(&fam_path, |line| {
        let c: Vec<&str> = crate::split_ws(line).collect();
        if c.len() >= 2 {
            ids.push(c[1].to_owned());
        }
        Ok(())
    })?;
    if ids.len() != people {
        return Err(Fail(format!(
            "{} lists {} people but ref.geno holds {people}: the model is damaged; run scripts/setup.sh --ancestry-fine-panel again",
            fam_path.display(),
            ids.len()
        )));
    }
    let samples_path = model
        .parent()
        .map_or_else(|| PathBuf::from("samples.tsv"), |p| p.join("samples.tsv"));
    let lines = crate::read_lines(&samples_path)?;
    let header: Vec<&str> = lines
        .first()
        .map(|l| l.split('\t').collect())
        .unwrap_or_default();
    if header != crate::panel::HEADER {
        return Err(Fail(format!(
            "{} does not start with the header of samples.tsv",
            samples_path.display()
        )));
    }
    let mut by_id: HashMap<&str, Vec<&str>> = HashMap::new();
    for line in &lines[1..] {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 4 {
            by_id.insert(f[0], f);
        }
    }
    for lm in levels.iter_mut() {
        lm.groups = vec![Vec::new(); lm.k];
        for (i, id) in ids.iter().enumerate() {
            let code = by_id.get(id.as_str()).map(|f| f[lm.level.column()]);
            let g = code.and_then(|c| lm.row.codes.iter().position(|x| x == c));
            let Some(g) = g else {
                return Err(Fail(format!(
                    "reference person {id} of ref.fam has no {} group of levels.tsv in {}: the model and the samples table do not match; run scripts/setup.sh --ancestry-fine-panel again",
                    lm.row.name,
                    samples_path.display()
                )));
            };
            lm.groups[g].push(i);
        }
    }
    Ok(())
}

/// The groups' A2 frequencies at the kept SNVs (`sites`, `ref.geno` rows)
/// with person i counted `people_w[i]` times; a SNV of weight 0 gets 0.5
/// in every group (it is never read).
fn group_freqs(
    geno: &RefGeno,
    sites: &[usize],
    snv_w: Option<&[u32]>,
    groups: &[Vec<usize>],
    people_w: Option<&[u32]>,
) -> Vec<f64> {
    let k = groups.len();
    // (person, weight, group) of everyone drawn at least once.
    let mut members: Vec<(usize, u32, usize)> = Vec::new();
    for (g, people) in groups.iter().enumerate() {
        for &i in people {
            let w = people_w.map_or(1, |w| w[i]);
            if w > 0 {
                members.push((i, w, g));
            }
        }
    }
    let mut out = vec![0.5; sites.len() * k];
    let (mut alt, mut called) = (vec![0u32; k], vec![0u32; k]);
    for (s, &site) in sites.iter().enumerate() {
        if snv_w.is_some_and(|w| w[s] == 0) {
            continue;
        }
        alt.iter_mut().for_each(|x| *x = 0);
        called.iter_mut().for_each(|x| *x = 0);
        for &(i, w, g) in &members {
            if let Some(c) = geno.get(site, i) {
                alt[g] += w * u32::from(c);
                called[g] += 2 * w;
            }
        }
        for g in 0..k {
            out[s * k + g] = a2_freq(alt[g], called[g]);
        }
    }
    out
}

/// What every bootstrap replicate reads.
struct Boot<'a> {
    seed: u64,
    blocks: &'a [Vec<usize>],
    genotypes: &'a [u8],
    levels: &'a [LevelModel],
    geno: &'a RefGeno,
    keep: &'a [usize],
    opts: &'a Options,
}

/// Replicate `r`'s fractions at every level.
fn replicate(b: &Boot<'_>, r: u64) -> Result<Vec<Vec<f64>>> {
    let mut rng = Rng::stream(b.seed, &format!("replicate {r}"));
    let snv_w = project::draw_weights(b.blocks, b.genotypes.len(), &mut rng);
    let mut out = Vec::with_capacity(b.levels.len());
    for lm in b.levels {
        let people_w = project::draw_people(&lm.groups, b.geno.people(), &mut rng);
        let freqs = group_freqs(b.geno, b.keep, Some(&snv_w), &lm.groups, Some(&people_w));
        let p = Problem::new(lm.k, b.genotypes, &freqs, Some(&snv_w)).map_err(|e| Fail(e.0))?;
        let fit = project::project(&p, b.opts).map_err(|e| {
            Fail(format!(
                "{} level, bootstrap replicate {r}: {}",
                lm.row.name, e.0
            ))
        })?;
        out.push(fit.q);
    }
    Ok(out)
}

fn q_line(q: &[f64]) -> String {
    let cells: Vec<String> = q.iter().map(|x| format!("{x:.6}")).collect();
    cells.join(" ") + "\n"
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let model = Path::new(a.req("--model"));
    let vcf = Path::new(a.req("--vcf"));
    let sample = a.req("--sample");
    if sample.is_empty() || sample.contains(['/', '\t', '\n']) {
        return Err(Fail(format!(
            "--sample {}: not a sample name",
            repr(sample)
        )));
    }
    let nonneg = |name: &str, default: i64| -> Result<i64> {
        let v = a.int(name, default).unwrap_or(default);
        if v < 0 {
            return Err(Fail(format!("{name} {v}: must be 0 or more")));
        }
        Ok(v)
    };
    let boot = nonneg("--bootstrap", 100)? as u64;
    let seed = nonneg("--seed", 1)? as u64;
    let threads = nonneg("--threads", 1)?.max(1) as usize;
    let min_sites = nonneg("--min-sites", 100)? as usize;
    let max_iterations = nonneg("--max-iterations", project::MAX_ITERATIONS as i64)? as usize;
    let opts = Options {
        max_iterations,
        ..Options::default()
    };

    let model_s = model.to_string_lossy().into_owned();
    let rows = read_levels(&model_s)?;
    let bim_path = model.join("ref.bim");
    let bim = read_bim(&bim_path)?;
    if bim.rows.is_empty() {
        return Err(Fail(format!("{} has no SNV", bim_path.display())));
    }
    let calls = genotypes_at(&bim, vcf)?;
    let (keep, genotypes): (Vec<usize>, Vec<u8>) = calls
        .iter()
        .enumerate()
        .filter_map(|(i, g)| g.map(|g| (i, g)))
        .unzip();
    if keep.len() < min_sites {
        return Err(Fail(format!(
            "only {} of the {} reference SNVs have a genotype in {} (at least {min_sites} needed)",
            keep.len(),
            bim.rows.len(),
            vcf.display()
        )));
    }
    let mut levels = Vec::new();
    for row in rows {
        let level = Level::parse(&row.name)
            .ok_or_else(|| Fail(format!("{model_s}/levels.tsv: unknown level {}", row.name)))?;
        let k = usize::try_from(row.k)
            .ok()
            .filter(|&k| k > 0 && k == row.codes.len())
            .ok_or_else(|| {
                Fail(format!(
                    "{model_s}/levels.tsv: level {} has K {} and {} group codes",
                    row.name,
                    row.k,
                    row.codes.len()
                ))
            })?;
        let freqs = read_p(
            &model.join(format!("{}.P", row.name)),
            k,
            &keep,
            bim.rows.len(),
        )?;
        levels.push(LevelModel {
            row,
            level,
            k,
            freqs,
            groups: Vec::new(),
        });
    }

    // The bootstrap's inputs: the blocks on the genetic map, and the
    // reference people of each group.
    let mut blocks = Vec::new();
    let mut geno = RefGeno::new(0);
    if boot > 0 {
        let geno_path = model.join("ref.geno");
        if !geno_path.is_file() {
            return Err(Fail(format!(
                "{} is missing: the model was made before og counted the group frequencies, so the bootstrap has no reference people to draw; run scripts/setup.sh --ancestry-fine-panel again",
                geno_path.display()
            )));
        }
        geno = RefGeno::load(&geno_path).map_err(|e| Fail(e.0))?;
        if geno.sites() != bim.rows.len() {
            return Err(Fail(format!(
                "{} holds {} SNVs but ref.bim {}: the model is damaged; run scripts/setup.sh --ancestry-fine-panel again",
                geno_path.display(),
                geno.sites(),
                bim.rows.len()
            )));
        }
        assign_groups(model, &mut levels, geno.people())?;
        // ref.geno must give the model's own frequencies, or the bootstrap
        // would measure the spread around another model.
        for lm in &levels {
            let full = group_freqs(&geno, &keep, None, &lm.groups, None);
            if let Some((i, d)) = full
                .iter()
                .zip(&lm.freqs)
                .map(|(x, y)| (x - y).abs())
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .filter(|&(_, d)| d > 1.000_001e-6)
            {
                return Err(Fail(format!(
                    "{} does not give {model_s}/{}.P (off by {d:.6} at SNV {} of ref.bim): the model is damaged; run scripts/setup.sh --ancestry-fine-panel again",
                    geno_path.display(),
                    lm.row.name,
                    keep[i / lm.k] + 1
                )));
            }
        }
        let maps = Path::new(a.get("--maps").ok_or_else(|| {
            Fail(
                "the bootstrap cuts blocks on the genetic map: give --maps (or --bootstrap 0)"
                    .into(),
            )
        })?);
        let mut names: HashMap<&str, (u32, GeneticMap)> = HashMap::new();
        let (mut chrom, mut cm) = (
            Vec::with_capacity(keep.len()),
            Vec::with_capacity(keep.len()),
        );
        for &i in &keep {
            let r = bim.rows[i];
            let c = bim.at(r.chrom);
            if !names.contains_key(c) {
                let path = maps.join(format!("chr{c}.map"));
                if !path.is_file() {
                    return Err(Fail(format!(
                        "{} is missing: the bootstrap cuts blocks on the genetic map; run scripts/setup.sh --ancestry-fine-panel again",
                        path.display()
                    )));
                }
                let map = GeneticMap::read(&path).map_err(|e| Fail(e.0))?;
                let id = u32::try_from(names.len()).expect("few chromosomes");
                names.insert(c, (id, map));
            }
            let (id, map) = &names[c];
            chrom.push(*id);
            cm.push(map.cm_at(r.pos));
        }
        blocks = project::blocks(&chrom, &cm, BLOCK_CM);
    }

    // The fractions, then the replicates on `threads` threads.
    let mut fits = Vec::with_capacity(levels.len());
    for lm in &levels {
        let p = Problem::new(lm.k, &genotypes, &lm.freqs, None).map_err(|e| Fail(e.0))?;
        let fit = project::project(&p, &opts)
            .map_err(|e| Fail(format!("{} level: {}", lm.row.name, e.0)))?;
        fits.push(fit);
    }
    let shared = Boot {
        seed,
        blocks: &blocks,
        genotypes: &genotypes,
        levels: &levels,
        geno: &geno,
        keep: &keep,
        opts: &opts,
    };
    let mut reps: Vec<Option<Result<Vec<Vec<f64>>>>> = (0..boot).map(|_| None).collect();
    std::thread::scope(|s| {
        let chunks: Vec<_> = (0..threads)
            .map(|t| {
                let shared = &shared;
                s.spawn(move || {
                    (1..=boot)
                        .filter(|r| (r - 1) % threads as u64 == t as u64)
                        .map(|r| (r, replicate(shared, r)))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        for h in chunks {
            for (r, res) in h.join().expect("a bootstrap thread panicked") {
                reps[(r - 1) as usize] = Some(res);
            }
        }
    });
    let reps: Vec<Vec<Vec<f64>>> = reps
        .into_iter()
        .map(|r| r.expect("every replicate ran"))
        .collect::<Result<_>>()?;

    let dir = PathBuf::from(a.req("--out-dir"));
    std::fs::create_dir_all(&dir).map_err(|e| crate::cannot_write(&dir, &e))?;
    let mut pending = Pending::new();
    let mut said = Vec::new();
    for (l, (lm, fit)) in levels.iter().zip(&fits).enumerate() {
        let base = format!("{sample}_{}", lm.row.name);
        pending.put(
            dir.join(format!("{base}.{}.Q", lm.k)),
            q_line(&fit.q).as_bytes(),
        )?;
        let boot_q: String = reps.iter().map(|r| q_line(&r[l])).collect();
        pending.put(dir.join(format!("{base}.boot.Q")), boot_q.as_bytes())?;
        let fit_tsv = format!(
            "key\tvalue\nlevel\t{}\ngroups\t{}\nsnvs\t{}\npasses\t{}\nloglik\t{:.6}\nkkt_slack\t{:e}\nkkt_gain\t{:e}\nblock_cm\t{BLOCK_CM}\nblocks\t{}\nreplicates\t{boot}\nreference_resampling\t{}\nseed\t{seed}\n",
            lm.row.name,
            lm.k,
            keep.len(),
            fit.iterations,
            fit.loglik,
            fit.slack,
            fit.gain,
            blocks.len(),
            if boot > 0 { "yes" } else { "no" }
        );
        pending.put(dir.join(format!("{base}.fit.tsv")), fit_tsv.as_bytes())?;
        said.push(format!(
            "{}: {} groups, {} passes, {boot} bootstrap replicates of {} blocks",
            lm.row.name,
            lm.k,
            fit.iterations,
            blocks.len()
        ));
    }
    let sites = format!(
        "{} of the {} SNVs of the model genotyped",
        keep.len(),
        bim.rows.len()
    );
    pending.put(
        dir.join(format!("{sample}_sites.txt")),
        format!("{sites}\n").as_bytes(),
    )?;
    pending.commit_all()?;
    let w = |e: std::io::Error| Fail(format!("cannot write stdout: {e}"));
    writeln!(out, "{sites}").map_err(w)?;
    for line in said {
        writeln!(out, "{line}").map_err(w)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    #[test]
    fn genotype_indexes_must_name_alleles() {
        assert_eq!(a2_copies(&["A", "G"], ["0", "1"], "G", "A"), Some(1));
        assert_eq!(a2_copies(&["A", "G"], ["0", "0"], "G", "A"), Some(2));
        assert_eq!(a2_copies(&["A", "G"], ["1", "1"], "G", "A"), Some(0));
        assert_eq!(a2_copies(&["A", "G"], ["0", "-1"], "G", "A"), None);
        assert_eq!(a2_copies(&["A"], ["0", "2"], "G", "A"), None);
        assert_eq!(a2_copies(&["A", "G"], ["x", "1"], "G", "A"), None);
        assert_eq!(a2_copies(&["A", "T"], ["0", "1"], "G", "A"), None);
        assert_eq!(gt_alleles("0/1/1"), None);
        assert_eq!(gt_alleles("./1"), None);
        assert_eq!(gt_alleles("1"), None);
    }

    /// A model of `snvs` SNVs on chr21 and chr22 (1 cM per Mb), K groups of
    /// `per` people each with frequencies 0.1 apart, and a person with
    /// fractions `q`.
    struct World {
        d: Dir,
    }

    const GROUPS: [&str; 3] = ["AFR", "WEA", "EAS"];

    fn world(name: &str, snvs: usize, per: usize, q: [f64; 3], seed: u64) -> World {
        let d = Dir::new(name);
        let mut rng = Rng::new(seed);
        let k = GROUPS.len();
        let people = k * per;
        let mut samples = String::from("sample\tpopulation\tcontinent\tregion\treference\tpaint\n");
        let mut fam = String::new();
        for group in GROUPS {
            for i in 0..per {
                let id = format!("{group}{i}");
                samples.push_str(&format!("{id}\tX\t{group}\t{group}\t1\t0\n"));
                fam.push_str(&format!("0\t{id}\t0\t0\t0\t-9\n"));
            }
        }
        let mut bim = String::new();
        let mut p_text = String::new();
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        );
        let mut geno = RefGeno::new(people);
        for j in 0..snvs {
            let chrom = if j < snvs / 2 { 21 } else { 22 };
            let pos = 1_000 + (j % (snvs / 2)) as i64 * 20_000;
            // The frequency of the ALT (A1) per group, then the people.
            let f: Vec<f64> = (0..k)
                .map(|_| {
                    let x = rng.uniform().powi(2);
                    if rng.below(2) == 0 { x } else { 1.0 - x }
                })
                .collect();
            let mut row = Vec::with_capacity(people);
            let mut counts = vec![(0u32, 0u32); k];
            for (g, fg) in f.iter().enumerate() {
                for _ in 0..per {
                    let c = u8::from(rng.uniform() < *fg) + u8::from(rng.uniform() < *fg);
                    counts[g].0 += u32::from(c);
                    counts[g].1 += 2;
                    row.push(Some(c));
                }
            }
            geno.push_site(&row).unwrap();
            let cells: Vec<String> = counts
                .iter()
                .map(|&(alt, called)| format!("{:.6}", a2_freq(alt, called)))
                .collect();
            p_text.push_str(&(cells.join(" ") + "\n"));
            bim.push_str(&format!("{chrom}\t{chrom}:{pos}\t0\t{pos}\tG\tA\n"));
            // The person's ALT copies from its fractions.
            let p_alt: f64 = q.iter().zip(&f).map(|(a, b)| a * b).sum();
            let alt = u8::from(rng.uniform() < p_alt) + u8::from(rng.uniform() < p_alt);
            let gt = ["0/0", "0/1", "1/1"][alt as usize];
            vcf.push_str(&format!("chr{chrom}\t{pos}\t.\tA\tG\t.\t.\t.\tGT\t{gt}\n"));
        }
        d.put("panel/samples.tsv", &samples);
        d.put("panel/model/ref.fam", &fam);
        d.put("panel/model/ref.bim", &bim);
        d.put("panel/model/continent.P", &p_text);
        d.put(
            "panel/model/levels.tsv",
            &format!("continent\t{k}\t{}\n", GROUPS.join(",")),
        );
        let mut bytes = Vec::new();
        geno.write(&mut bytes).unwrap();
        std::fs::write(d.0.join("panel/model/ref.geno"), bytes).unwrap();
        for c in [21, 22] {
            d.put(
                &format!("panel/maps/chr{c}.map"),
                &format!("chr{c} . 0.0 1\nchr{c} . 100.0 100000001\n"),
            );
        }
        d.put("g.vcf", &vcf);
        World { d }
    }

    impl World {
        fn project(&self, extra: &[&str]) -> (u8, String, String) {
            let (model, vcf, maps, out) = (
                self.d.path("panel/model"),
                self.d.path("g.vcf"),
                self.d.path("panel/maps"),
                self.d.path("out"),
            );
            let mut args = vec![
                "project",
                "--model",
                &model,
                "--vcf",
                &vcf,
                "--sample",
                "S",
                "--maps",
                &maps,
                "--out-dir",
                &out,
            ];
            args.extend_from_slice(extra);
            run(&args)
        }

        fn q(&self, name: &str) -> Vec<Vec<f64>> {
            self.d
                .read(&format!("out/{name}"))
                .lines()
                .map(|l| l.split(' ').map(|x| x.parse().unwrap()).collect())
                .collect()
        }
    }

    fn spread(rows: &[Vec<f64>], k: usize) -> f64 {
        let col: Vec<f64> = rows.iter().map(|r| r[k]).collect();
        let lo = crate::report::percentile(&col, 0.025).unwrap();
        let hi = crate::report::percentile(&col, 0.975).unwrap();
        hi - lo
    }

    #[test]
    fn the_command_projects_and_bootstraps() {
        let w = world("proj-cmd", 2_000, 30, [0.7, 0.2, 0.1], 5);
        let (rc, out, err) = w.project(&["--bootstrap", "20", "--threads", "3"]);
        assert_eq!((rc, err.as_str()), (0, ""));
        assert!(
            out.starts_with("2000 of the 2000 SNVs of the model genotyped\ncontinent: 3 groups, "),
            "{out}"
        );
        let q = &w.q("S_continent.3.Q")[0];
        for (got, want) in q.iter().zip([0.7, 0.2, 0.1]) {
            assert!((got - want).abs() < 0.05, "{q:?}");
        }
        let boot = w.q("S_continent.boot.Q");
        assert_eq!(boot.len(), 20);
        assert!(
            boot.iter()
                .all(|r| r.len() == 3 && (r.iter().sum::<f64>() - 1.0).abs() < 1e-5)
        );
        let fit = w.d.read("out/S_continent.fit.tsv");
        // 2 chromosomes of 20 Mb at 1 cM per Mb: 4 blocks of 5 cM each, 8 in all.
        assert!(
            fit.contains("\nblocks\t8\n") && fit.contains("\nreplicates\t20\n"),
            "{fit}"
        );
        assert_eq!(
            w.d.read("out/S_sites.txt"),
            "2000 of the 2000 SNVs of the model genotyped\n"
        );
        // The same seed gives the same bytes, whatever the threads.
        let first = w.d.read("out/S_continent.boot.Q");
        let (rc, _, _) = w.project(&["--bootstrap", "20", "--threads", "1"]);
        assert_eq!(rc, 0);
        let again = w.d.read("out/S_continent.boot.Q");
        assert_eq!(first, again);
        // And another seed gives other replicates.
        let (rc, _, _) = w.project(&["--bootstrap", "20", "--seed", "2"]);
        assert_eq!(rc, 0);
        assert_ne!(first, w.d.read("out/S_continent.boot.Q"));
    }

    #[test]
    fn drawing_the_people_again_widens_a_small_group() {
        // A group of 10 people: its frequencies are uncertain, and the
        // range with the people drawn again is wider than with the model's
        // frequencies fixed (the same blocks, the same seed).
        let w = world("proj-small", 2_000, 10, [0.5, 0.3, 0.2], 9);
        let (rc, _, err) = w.project(&["--bootstrap", "60"]);
        assert_eq!((rc, err.as_str()), (0, ""));
        let resampled = w.q("S_continent.boot.Q");
        // Fixed frequencies: the same block draws, through the engine.
        let model = w.d.0.join("panel/model");
        let bim = read_bim(&model.join("ref.bim")).unwrap();
        let calls = genotypes_at(&bim, &w.d.0.join("g.vcf")).unwrap();
        let (keep, g): (Vec<usize>, Vec<u8>) = calls
            .iter()
            .enumerate()
            .filter_map(|(i, c)| c.map(|c| (i, c)))
            .unzip();
        let f = read_p(&model.join("continent.P"), 3, &keep, bim.rows.len()).unwrap();
        let chrom: Vec<u32> = keep.iter().map(|&i| u32::from(i >= 1_000)).collect();
        let cm: Vec<f64> = keep.iter().map(|&i| bim.rows[i].pos as f64 / 1e6).collect();
        let blocks = project::blocks(&chrom, &cm, BLOCK_CM);
        let fixed: Vec<Vec<f64>> = (1..=60u64)
            .map(|r| {
                let mut rng = Rng::stream(1, &format!("replicate {r}"));
                let sw = project::draw_weights(&blocks, g.len(), &mut rng);
                let p = Problem::new(3, &g, &f, Some(&sw)).unwrap();
                project::project(&p, &Options::default()).unwrap().q
            })
            .collect();
        let (wide, narrow) = (spread(&resampled, 0), spread(&fixed, 0));
        assert!(
            wide > narrow * 1.2,
            "with people drawn again {wide}, fixed {narrow}"
        );
    }

    #[test]
    fn what_stops_the_command() {
        let w = world("proj-stop", 400, 10, [0.6, 0.3, 0.1], 3);
        let (rc, _, err) = w.project(&["--min-sites", "500", "--bootstrap", "0"]);
        assert_eq!(rc, 1);
        assert!(
            err.contains("only 400 of the 400 reference SNVs have a genotype in"),
            "{err}"
        );
        // Without ref.geno: no bootstrap, one line to run setup again; with
        // --bootstrap 0 the fractions still come.
        std::fs::remove_file(w.d.0.join("panel/model/ref.geno")).unwrap();
        let (rc, _, err) = w.project(&["--bootstrap", "5"]);
        assert_eq!(rc, 1);
        assert!(
            err.starts_with("ERROR: ")
                && err.ends_with("run scripts/setup.sh --ancestry-fine-panel again\n")
                && err.lines().count() == 1,
            "{err}"
        );
        assert!(!w.d.0.join("out").exists(), "a failure writes nothing");
        let (rc, _, err) = w.project(&["--bootstrap", "0"]);
        assert_eq!((rc, err.as_str()), (0, ""));
        assert_eq!(w.d.read("out/S_continent.boot.Q"), "");
        // A .P of the wrong length.
        w.d.put("panel/model/continent.P", "0.5 0.5 0.5\n");
        let (rc, _, err) = w.project(&["--bootstrap", "0"]);
        assert_eq!(rc, 1);
        assert!(
            err.contains("continent.P has 1 rows but ref.bim 400 SNVs"),
            "{err}"
        );
        // A truncated run fails its KKT check, and says so.
        let w = world("proj-kkt", 400, 10, [0.6, 0.3, 0.1], 3);
        let (rc, _, err) = w.project(&["--bootstrap", "0", "--max-iterations", "3"]);
        assert_eq!(rc, 1);
        assert!(
            err.contains("continent level: the projection stopped after 3 passes"),
            "{err}"
        );
    }

    #[test]
    fn a_damaged_model_is_refused() {
        // ref.geno that does not give the .P: one person's genotypes
        // changed.
        let w = world("proj-damaged", 400, 10, [0.6, 0.3, 0.1], 4);
        let mut bytes = std::fs::read(w.d.0.join("panel/model/ref.geno")).unwrap();
        let at = 28 + 5 * 8; // SNV 6, people 1 to 4 (8 bytes a row for 30 people)
        bytes[at] ^= 0b0101_0101;
        std::fs::write(w.d.0.join("panel/model/ref.geno"), bytes).unwrap();
        let (rc, _, err) = w.project(&["--bootstrap", "2"]);
        assert_eq!(rc, 1);
        assert!(
            err.contains("does not give") && err.contains("the model is damaged"),
            "{err}"
        );
    }
}
