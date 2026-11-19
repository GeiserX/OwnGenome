//! Local-ancestry painting, our own engine in place of RFMix: a Li and
//! Stephens copying model whose hidden state is (ancestry, reference
//! haplotype), as Browning, Waples and Browning published it for FLARE
//! (AJHG 2023). The code is written from the published equations; no FLARE
//! or RFMix source was read.
//!
//! Along a target haplotype, ancestry events come at rate T per Morgan, and
//! each draws the new ancestry from the proportions mu (it may be the old
//! one) and a reference haplotype of it uniformly. Inside an ancestry a the
//! copied haplotype jumps at rate rho_a per Morgan to one drawn uniformly
//! from the n_a haplotypes of a. The copied haplotype carries the target's
//! allele with chance 1 - epsilon. [`hmm`] has the equations; each forward
//! and backward step costs O(H) through per-group totals, and the forward
//! values are kept only every sqrt(M) sites.
//!
//! [`paint`] fits mu, T, epsilon and rho_a by expectation-maximisation over
//! the whole genome, every chromosome and haplotype at once, for up to
//! [`Options::rounds`] rounds, and calls each site's ancestry as the group
//! of highest posterior. The start is [`Params::initial`]: mu from the
//! continental projection, T = 10, epsilon = 0.01 and rho_a = 4 Ne / n_a.
//! Chromosomes and haplotypes run on [`Options::threads`] threads; the
//! sums are taken in a fixed order, so any thread count gives the same
//! bytes.
//!
//! Memory: the reference haplotypes take M H / 8 bytes over all M sites of
//! the genome, each running haplotype about 16 sqrt(M_c) H bytes for its
//! chromosome of M_c sites, and the posteriors 8 A bytes per site and
//! haplotype for A groups.

pub mod bits;
pub mod hmm;

pub use bits::{Bits, RefPanel};
pub use hmm::{Fb, Stats, checkpoint_block, forward_backward};

use crate::simulate::MAX_PER_MORGAN;
use crate::{Error, Result, repr};
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The effective population size the start of rho_a takes, 4 Ne / n_a per
/// Morgan.
pub const DEFAULT_NE: f64 = 10_000.0;
/// The start of T, ancestry events per Morgan.
pub const START_T: f64 = 10.0;
/// The start of epsilon.
pub const START_EPSILON: f64 = 0.01;
/// The most EM rounds [`paint`] runs by default.
pub const DEFAULT_ROUNDS: usize = 5;
/// EM stops early when no parameter moves by more than this share.
pub const TOLERANCE: f64 = 1e-4;

/// The model's parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct Params {
    /// The proportion of each group, summing to 1.
    pub mu: Vec<f64>,
    /// T, ancestry events per Morgan.
    pub generations: f64,
    /// The chance that the copied haplotype does not carry the allele.
    pub epsilon: f64,
    /// rho_a, copying jumps per Morgan inside each group.
    pub rho: Vec<f64>,
}

impl Params {
    /// The start: `mu` (the continental projection), T = 10, epsilon =
    /// 0.01 and rho_a = 4 `ne` / n_a, `sizes` holding n_a.
    pub fn initial(mu: Vec<f64>, sizes: &[usize], ne: f64) -> Params {
        Params {
            mu,
            generations: START_T,
            epsilon: START_EPSILON,
            rho: sizes.iter().map(|&n| 4.0 * ne / n as f64).collect(),
        }
    }

    /// Checks the numbers for `groups` groups: mu above 0 and summing to 1
    /// (within 1e-6), T above 0 and at most [`MAX_PER_MORGAN`], epsilon
    /// above 0 and at most 0.5, every rho_a finite and above 0.
    pub fn check(&self, groups: usize) -> Result<()> {
        if self.mu.len() != groups || self.rho.len() != groups {
            return Err(Error(format!(
                "{} proportions and {} jump rates for {groups} groups",
                self.mu.len(),
                self.rho.len()
            )));
        }
        if let Some(w) = self.mu.iter().find(|w| !(w.is_finite() && **w > 0.0)) {
            return Err(Error(format!("a proportion is {w}: each must be above 0")));
        }
        let total: f64 = self.mu.iter().sum();
        if (total - 1.0).abs() > 1e-6 {
            return Err(Error(format!(
                "the ancestry proportions sum to {total}, not 1"
            )));
        }
        let t = self.generations;
        if !(t.is_finite() && t > 0.0 && t <= MAX_PER_MORGAN) {
            return Err(Error(format!(
                "T is {t}: ancestry events per Morgan must be above 0 and at most {MAX_PER_MORGAN}"
            )));
        }
        let eps = self.epsilon;
        if !(eps.is_finite() && eps > 0.0 && eps <= 0.5) {
            return Err(Error(format!(
                "epsilon is {eps}: it must be above 0 and at most 0.5"
            )));
        }
        if let Some(r) = self.rho.iter().find(|r| !(r.is_finite() && **r > 0.0)) {
            return Err(Error(format!(
                "a copying jump rate is {r}: each must be above 0"
            )));
        }
        Ok(())
    }

    /// The next parameters from the expected counts of a pass.
    fn update(&self, s: &Stats) -> Params {
        let groups = self.mu.len();
        // mu: the mean posterior of each ancestry, kept above 1e-6 so no
        // group drops out of the model.
        let mut mu: Vec<f64> = s.ancestry.iter().map(|x| (x / s.sites).max(1e-6)).collect();
        let total: f64 = mu.iter().sum();
        mu.iter_mut().for_each(|x| *x /= total);
        // T: expected events per Morgan walked.
        let generations = if s.morgans > 0.0 {
            (s.events / s.morgans).clamp(1e-3, MAX_PER_MORGAN)
        } else {
            self.generations
        };
        let epsilon = (s.mismatches / s.sites).clamp(1e-6, 0.5);
        let rho = (0..groups)
            .map(|a| {
                if s.exposure[a] > 0.0 {
                    (s.jumps[a] / s.exposure[a]).clamp(1e-3, 1e7)
                } else {
                    self.rho[a]
                }
            })
            .collect();
        Params {
            mu,
            generations,
            epsilon,
            rho,
        }
    }

    /// The largest share by which any parameter differs from `other`;
    /// infinite when a share is not a number, so a NaN never reads as
    /// converged.
    fn moved(&self, other: &Params) -> f64 {
        let rel = |a: f64, b: f64| {
            let r = (a - b).abs() / b.abs().max(f64::MIN_POSITIVE);
            if r.is_nan() { f64::INFINITY } else { r }
        };
        let mut m = rel(self.generations, other.generations).max(rel(self.epsilon, other.epsilon));
        for (a, b) in self
            .mu
            .iter()
            .zip(&other.mu)
            .chain(self.rho.iter().zip(&other.rho))
        {
            m = m.max(rel(*a, *b));
        }
        m
    }
}

/// One chromosome to paint.
#[derive(Debug, Clone, PartialEq)]
pub struct ChromInput {
    pub name: String,
    /// The genetic position of each site, cM, never going down.
    pub cm: Vec<f64>,
    pub panel: RefPanel,
    /// The target's haplotypes at the panel's sites.
    pub targets: Vec<Bits>,
}

/// How [`paint`] runs.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Worker threads, at least 1.
    pub threads: usize,
    /// The most EM rounds; 0 paints with the start parameters.
    pub rounds: usize,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            threads: 1,
            rounds: DEFAULT_ROUNDS,
        }
    }
}

/// One painted haplotype.
#[derive(Debug, Clone, PartialEq)]
pub struct PaintedHap {
    /// Site by site, the posterior of each group: `posterior[m * A + a]`.
    pub posterior: Vec<f64>,
    /// Site by site, the group of highest posterior (the first on a tie).
    pub calls: Vec<usize>,
}

/// One painted chromosome.
#[derive(Debug, Clone, PartialEq)]
pub struct PaintedChrom {
    pub name: String,
    pub haplotypes: Vec<PaintedHap>,
}

/// The painting of a genome.
#[derive(Debug, Clone, PartialEq)]
pub struct Painting {
    /// The parameters the painting was made with.
    pub params: Params,
    /// The EM rounds run (parameter updates made).
    pub rounds: usize,
    /// The log-likelihood of the genome under [`Painting::params`].
    pub loglik: f64,
    pub chroms: Vec<PaintedChrom>,
}

fn check_input(chroms: &[ChromInput], start: &Params, opts: &Options) -> Result<()> {
    if opts.threads == 0 {
        return Err(Error("the painting needs at least one thread".into()));
    }
    let Some(first) = chroms.first() else {
        return Err(Error("no chromosome to paint".into()));
    };
    start.check(first.panel.groups())?;
    if chroms.iter().all(|c| c.targets.is_empty()) {
        return Err(Error("no target haplotype to paint".into()));
    }
    for c in chroms {
        let name = repr(&c.name);
        if c.panel.sizes() != first.panel.sizes() {
            return Err(Error(format!(
                "chromosome {name} has another reference panel than {}",
                repr(&first.name)
            )));
        }
        let m = c.panel.sites();
        if m == 0 {
            return Err(Error(format!("chromosome {name} has no site")));
        }
        if c.cm.len() != m {
            return Err(Error(format!(
                "chromosome {name}: {} map positions for {m} sites",
                c.cm.len()
            )));
        }
        if let Some(i) = (0..m).find(|&i| !c.cm[i].is_finite() || (i > 0 && c.cm[i] < c.cm[i - 1]))
        {
            return Err(Error(format!(
                "chromosome {name}: the map position of site {i} is {} (it must be finite and never go down)",
                c.cm[i]
            )));
        }
        if let Some(t) = c.targets.iter().position(|t| t.len() != m) {
            return Err(Error(format!(
                "chromosome {name}: target haplotype {t} has {} sites, the panel {m}",
                c.targets[t].len()
            )));
        }
    }
    Ok(())
}

/// Runs one forward-backward per (chromosome, haplotype) on `threads`
/// threads; the results in that order.
fn pass(chroms: &[ChromInput], params: &Params, threads: usize) -> Vec<Fb> {
    let units: Vec<(usize, usize)> = chroms
        .iter()
        .enumerate()
        .flat_map(|(c, ch)| (0..ch.targets.len()).map(move |t| (c, t)))
        .collect();
    let run = |&(c, t): &(usize, usize)| {
        let ch = &chroms[c];
        forward_backward(
            &ch.panel,
            &ch.cm,
            &ch.targets[t],
            params,
            checkpoint_block(ch.panel.sites()),
        )
    };
    if threads <= 1 || units.len() <= 1 {
        return units.iter().map(run).collect();
    }
    let next = AtomicUsize::new(0);
    let mut done: Vec<(usize, Fb)> = std::thread::scope(|s| {
        let workers: Vec<_> = (0..threads.min(units.len()))
            .map(|_| {
                s.spawn(|| {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(u) = units.get(i) else { break };
                        out.push((i, run(u)));
                    }
                    out
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().expect("a painting thread panicked"))
            .collect()
    });
    done.sort_by_key(|(i, _)| *i);
    done.into_iter().map(|(_, fb)| fb).collect()
}

/// Paints every target haplotype of `chroms`, fitting the parameters from
/// `start` by EM for up to `opts.rounds` rounds over the whole genome. A
/// round stops the fit early when no parameter moves by more than
/// [`TOLERANCE`]; the painting is then made with the parameters before
/// that last update.
pub fn paint(chroms: &[ChromInput], start: Params, opts: &Options) -> Result<Painting> {
    check_input(chroms, &start, opts)?;
    let groups = chroms[0].panel.groups();
    let mut params = start;
    let mut rounds = 0;
    let (fbs, total) = loop {
        let fbs = pass(chroms, &params, opts.threads);
        let mut total = Stats::zero(groups);
        for fb in &fbs {
            total.add(&fb.stats);
        }
        if rounds == opts.rounds {
            break (fbs, total);
        }
        let next = params.update(&total);
        next.check(groups)?;
        if next.moved(&params) <= TOLERANCE {
            break (fbs, total);
        }
        params = next;
        rounds += 1;
    };
    let mut fbs = fbs.into_iter();
    let out = chroms
        .iter()
        .map(|c| PaintedChrom {
            name: c.name.clone(),
            haplotypes: (0..c.targets.len())
                .map(|_| {
                    let fb = fbs.next().expect("one pass per haplotype");
                    let calls = fb
                        .posterior
                        .chunks(groups)
                        .map(|p| {
                            (1..groups).fold(0, |best, a| if p[a] > p[best] { a } else { best })
                        })
                        .collect();
                    PaintedHap {
                        posterior: fb.posterior,
                        calls,
                    }
                })
                .collect(),
        })
        .collect();
    Ok(Painting {
        params,
        rounds,
        loglik: total.loglik,
        chroms: out,
    })
}

/// Writes `paint/params.tsv`: a header, then one `parameter group value`
/// row per value (`.` for a parameter of no group): T, epsilon, mu and rho
/// of each group, the EM rounds and the log-likelihood. `names` names the
/// groups in their order.
pub fn write_params(out: &mut dyn Write, painting: &Painting, names: &[String]) -> Result<()> {
    let p = &painting.params;
    if names.len() != p.mu.len() {
        return Err(Error(format!(
            "{} group names for {} groups",
            names.len(),
            p.mu.len()
        )));
    }
    let mut text = String::from("#parameter\tgroup\tvalue\n");
    text.push_str(&format!("T\t.\t{}\n", p.generations));
    text.push_str(&format!("epsilon\t.\t{}\n", p.epsilon));
    for (g, w) in names.iter().zip(&p.mu) {
        text.push_str(&format!("mu\t{g}\t{w}\n"));
    }
    for (g, r) in names.iter().zip(&p.rho) {
        text.push_str(&format!("rho\t{g}\t{r}\n"));
    }
    text.push_str(&format!("em_rounds\t.\t{}\n", painting.rounds));
    text.push_str(&format!("loglik\t.\t{}\n", painting.loglik));
    out.write_all(text.as_bytes())
        .map_err(|e| Error(format!("cannot write the painting parameters: {e}")))
}
