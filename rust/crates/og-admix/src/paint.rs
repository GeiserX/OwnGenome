//! Local-ancestry painting, our own engine in place of RFMix, on two
//! levels.
//!
//! Level 1, per group a ([`split::local_evidence`]): a Li and Stephens
//! copying model over the n_a reference haplotypes of a alone, as
//! Browning, Waples and Browning published it for FLARE (AJHG 2023): the
//! copied haplotype jumps at rate rho_a per Morgan to one drawn uniformly
//! from a, and carries the target's allele with chance 1 - epsilon. Run
//! left to right and right to left, its per-site predictive probabilities
//! give the local evidence ell_a(m) of each group. Each group's copying
//! keeps running whatever the ancestry is, so a switch of copied haplotype
//! inside a group never costs an ancestry switch.
//!
//! Level 2 ([`split::ancestry_fb`]): an ancestry HMM over the A groups.
//! Ancestry events come at rate T per Morgan and draw the new ancestry
//! from the proportions mu (it may be the old one); the emission at site m
//! is exp(lambda ell_a(m)), lambda = [`split::LAMBDA`]. lambda below 1
//! tempers the evidence, whose per-site terms inside one copied segment
//! are strongly correlated.
//!
//! [`paint`] computes level 1 once, with epsilon = 0.01 and rho_a = 4 Ne /
//! n_a, then fits T alone by expectation-maximisation over the whole
//! genome for up to [`Options::rounds`] rounds. mu stays at the start, the
//! continental projection: fitted from the posterior it drifts towards
//! groups that are not there. Each site's ancestry is the group of highest
//! posterior. Chromosomes and haplotypes run on [`Options::threads`]
//! threads; the sums are taken in a fixed order, so any thread count gives
//! the same bytes.
//!
//! [`hmm`] keeps the joint (ancestry, reference haplotype) model of FLARE,
//! which over-segments: one sharp emission picks both the copied haplotype
//! and the ancestry, so a few sites where another group's haplotype matches
//! better buy an ancestry switch, and EM then raises T to explain them.
//!
//! Memory: the reference haplotypes take M H / 8 bytes over all M sites of
//! the genome, the level-1 evidence 8 A bytes per site and haplotype, and
//! each running haplotype 24 A bytes per site of its chromosome.

pub mod bits;
pub mod hmm;
pub mod split;

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

    /// The next parameters from the expected counts of a pass: T, expected
    /// events per Morgan walked. mu, epsilon and rho_a stay.
    fn update(&self, s: &Stats) -> Params {
        let generations = if s.morgans > 0.0 {
            (s.events / s.morgans).clamp(1e-3, MAX_PER_MORGAN)
        } else {
            self.generations
        };
        Params {
            generations,
            ..self.clone()
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
    /// The log-likelihood of the genome under the ancestry level with
    /// [`Painting::params`], its evidence tempered by [`split::LAMBDA`]: a
    /// pseudo-likelihood, comparable between runs of this engine only.
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

/// Runs `f` on every item on `threads` threads; the results in order.
fn par_map<T: Sync, R: Send>(items: &[T], threads: usize, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    if threads <= 1 || items.len() <= 1 {
        return items.iter().map(&f).collect();
    }
    let next = AtomicUsize::new(0);
    let mut done: Vec<(usize, R)> = std::thread::scope(|s| {
        let workers: Vec<_> = (0..threads.min(items.len()))
            .map(|_| {
                s.spawn(|| {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(u) = items.get(i) else { break };
                        out.push((i, f(u)));
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
    done.into_iter().map(|(_, r)| r).collect()
}

/// Paints every target haplotype of `chroms`, fitting the parameters from
/// `start` by EM for up to `opts.rounds` rounds over the whole genome. A
/// round stops the fit early when no parameter moves by more than
/// [`TOLERANCE`]; the painting is then made with the parameters before
/// that last update.
pub fn paint(chroms: &[ChromInput], start: Params, opts: &Options) -> Result<Painting> {
    check_input(chroms, &start, opts)?;
    let groups = chroms[0].panel.groups();
    // One unit per (chromosome, haplotype), in that order.
    let units: Vec<(usize, usize)> = chroms
        .iter()
        .enumerate()
        .flat_map(|(c, ch)| (0..ch.targets.len()).map(move |t| (c, t)))
        .collect();
    // Level 1 depends on epsilon and rho_a only, which the fit keeps.
    let evidence: Vec<Vec<f64>> = par_map(&units, opts.threads, |&(c, t)| {
        let ch = &chroms[c];
        split::local_evidence(&ch.panel, &ch.cm, &ch.targets[t], &start)
    });
    let index: Vec<usize> = (0..units.len()).collect();
    let mut params = start;
    let mut rounds = 0;
    let (fbs, total) = loop {
        let fbs = par_map(&index, opts.threads, |&i| {
            split::ancestry_fb(
                &evidence[i],
                &chroms[units[i].0].cm,
                &params,
                split::LAMBDA,
                groups,
            )
        });
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
/// of each group, the EM rounds and the log-likelihood ([`Painting::loglik`]).
/// `names` names the groups in their order.
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
