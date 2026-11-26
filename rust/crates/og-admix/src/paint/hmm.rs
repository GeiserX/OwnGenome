//! Forward-backward over the (ancestry, reference haplotype) states of one
//! target haplotype on one chromosome, in O(H) per site, with the forward
//! values kept only at checkpoints.
//!
//! Between sites m-1 and m, d Morgans apart, with pT = 1 - exp(-T d) and
//! pR(a) = 1 - exp(-rho_a d), the chance of moving from state (a, h) to
//! state (b, k) is
//!
//! ```text
//! P = (1 - pT) [ (1 - pR(a)) [h = k] + pR(a) / n_a ] [a = b]
//!   + pT mu_b / n_b
//! ```
//!
//! The first term keeps the ancestry (and, unless the copying jumps, the
//! haplotype); the second is an ancestry event, whose new ancestry is drawn
//! from mu (it may be the old one) and whose haplotype is drawn uniformly in
//! it. Every term is either the state itself or a group total, so with the
//! per-group sums S_a = sum of alpha(a, h) over h the forward step is
//!
//! ```text
//! alpha_m(b, k) = e_m(k) [ (1 - pT)(1 - pR(b)) alpha_{m-1}(b, k)
//!                        + (1 - pT) pR(b) S_b / n_b + pT mu_b / n_b sum(S) ]
//! ```
//!
//! and the backward step, with g(b, k) = e_m(k) beta_m(b, k) and its group
//! sums G_b, is
//!
//! ```text
//! beta_{m-1}(a, h) = (1 - pT)(1 - pR(a)) g(a, h) + (1 - pT) pR(a) G_a / n_a
//!                  + pT sum_b mu_b G_b / n_b
//! ```
//!
//! The emission e_m(k) is 1 - epsilon when reference haplotype k carries
//! the target's allele at site m and epsilon when it does not. Both passes
//! rescale each site to sum 1; the forward scales give the log-likelihood.
//!
//! The forward values of every `block`-th site are kept. The backward pass
//! recomputes one block at a time from its checkpoint, with the same code
//! and in the same order, so the values are the same bits as in the first
//! pass. With `block` near sqrt(M) the memory is about 2 sqrt(M) H values
//! instead of M H.

use super::Params;
use super::bits::{Bits, RefPanel};

/// What one haplotype's pass adds to the EM step: expected counts, summed
/// over the sites of one chromosome.
#[derive(Debug, Clone, PartialEq)]
pub struct Stats {
    /// Sites.
    pub sites: f64,
    /// The expected number of sites where the copied haplotype does not
    /// carry the target's allele.
    pub mismatches: f64,
    /// Per group, the summed posterior of that ancestry over the sites.
    pub ancestry: Vec<f64>,
    /// The expected number of ancestry events.
    pub events: f64,
    /// Morgans from the first site to the last.
    pub morgans: f64,
    /// Per group, the expected number of copying jumps inside that
    /// ancestry.
    pub jumps: Vec<f64>,
    /// Per group, the expected Morgans spent in that ancestry (each
    /// interval weighted by the ancestry's posterior at its left end).
    pub exposure: Vec<f64>,
    /// The log-likelihood of the haplotype.
    pub loglik: f64,
}

impl Stats {
    pub fn zero(groups: usize) -> Stats {
        Stats {
            sites: 0.0,
            mismatches: 0.0,
            ancestry: vec![0.0; groups],
            events: 0.0,
            morgans: 0.0,
            jumps: vec![0.0; groups],
            exposure: vec![0.0; groups],
            loglik: 0.0,
        }
    }

    /// Adds `other` in place.
    pub fn add(&mut self, other: &Stats) {
        self.sites += other.sites;
        self.mismatches += other.mismatches;
        self.events += other.events;
        self.morgans += other.morgans;
        self.loglik += other.loglik;
        for (a, b) in [
            (&mut self.ancestry, &other.ancestry),
            (&mut self.jumps, &other.jumps),
            (&mut self.exposure, &other.exposure),
        ] {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
    }
}

/// The result of one haplotype's pass.
#[derive(Debug, Clone, PartialEq)]
pub struct Fb {
    /// Site by site, the posterior of each group: `posterior[m * A + a]`.
    pub posterior: Vec<f64>,
    pub stats: Stats,
}

/// The transition terms of one interval, per group.
struct Step {
    p_t: f64,
    /// (1 - pT)(1 - pR(a)).
    stay: Vec<f64>,
    /// (1 - pT) pR(a) / n_a.
    jump: Vec<f64>,
    /// pT mu_a / n_a.
    enter: Vec<f64>,
}

impl Step {
    fn new(groups: usize) -> Step {
        Step {
            p_t: 0.0,
            stay: vec![0.0; groups],
            jump: vec![0.0; groups],
            enter: vec![0.0; groups],
        }
    }

    fn set(&mut self, params: &Params, sizes: &[usize], morgans: f64) {
        let p_t = -(-params.generations * morgans).exp_m1();
        self.p_t = p_t;
        for (a, &n) in sizes.iter().enumerate() {
            let n = n as f64;
            let p_r = -(-params.rho[a] * morgans).exp_m1();
            self.stay[a] = (1.0 - p_t) * (1.0 - p_r);
            self.jump[a] = (1.0 - p_t) * p_r / n;
            self.enter[a] = p_t * params.mu[a] / n;
        }
    }
}

struct Engine<'a> {
    panel: &'a RefPanel,
    cm: &'a [f64],
    target: &'a Bits,
    params: &'a Params,
    step: Step,
    sums: Vec<f64>,
}

impl Engine<'_> {
    fn emission(&self, m: usize, h: usize) -> f64 {
        if self.panel.allele(m, h) == self.target.get(m) {
            1.0 - self.params.epsilon
        } else {
            self.params.epsilon
        }
    }

    fn morgans(&self, m: usize) -> f64 {
        (self.cm[m] - self.cm[m - 1]) / 100.0
    }

    /// Writes the forward values of site `m` into `out`, scaled to sum 1,
    /// from those of site m-1 (`prev`, scaled to sum 1); the log of the
    /// scale.
    fn forward(&mut self, m: usize, prev: Option<&[f64]>, out: &mut [f64]) -> f64 {
        let group = self.panel.group_all();
        let sizes = self.panel.sizes();
        match prev {
            None => {
                for (h, o) in out.iter_mut().enumerate() {
                    let g = group[h];
                    *o = self.emission(0, h) * self.params.mu[g] / sizes[g] as f64;
                }
            }
            Some(prev) => {
                let d = self.morgans(m);
                self.step.set(self.params, sizes, d);
                self.sums.iter_mut().for_each(|s| *s = 0.0);
                for (h, &p) in prev.iter().enumerate() {
                    self.sums[group[h]] += p;
                }
                let total: f64 = self.sums.iter().sum();
                for (a, s) in self.sums.iter_mut().enumerate() {
                    // Everything a state of group a receives that does not
                    // depend on its own previous value.
                    *s = self.step.jump[a] * *s + self.step.enter[a] * total;
                }
                for (h, o) in out.iter_mut().enumerate() {
                    let g = group[h];
                    *o = self.emission(m, h) * (self.step.stay[g] * prev[h] + self.sums[g]);
                }
            }
        }
        let c: f64 = out.iter().sum();
        out.iter_mut().for_each(|o| *o /= c);
        c.ln()
    }
}

/// The checkpoint spacing for `sites` sites: ceil(sqrt(sites)), at least 1.
pub fn checkpoint_block(sites: usize) -> usize {
    sites.isqrt().max(1) + usize::from(sites.isqrt().pow(2) < sites)
}

/// Runs forward-backward for `target` against `panel`, `cm` the genetic
/// position of each site, keeping the forward values every `block` sites
/// (`block` at least the number of sites keeps them all).
pub fn forward_backward(
    panel: &RefPanel,
    cm: &[f64],
    target: &Bits,
    params: &Params,
    block: usize,
) -> Fb {
    let sites = panel.sites();
    let hs = panel.haplotypes();
    let groups = panel.groups();
    assert!(sites > 0 && cm.len() == sites && target.len() == sites);
    let k = block.clamp(1, sites);
    let blocks = sites.div_ceil(k);
    let mut e = Engine {
        panel,
        cm,
        target,
        params,
        step: Step::new(groups),
        sums: vec![0.0; groups],
    };
    let mut stats = Stats::zero(groups);

    // Forward, keeping every k-th site.
    let mut checkpoints = vec![0.0; blocks * hs];
    let (mut prev, mut cur) = (vec![0.0; hs], vec![0.0; hs]);
    for m in 0..sites {
        stats.loglik += e.forward(m, (m > 0).then_some(&prev[..]), &mut cur);
        if m % k == 0 {
            checkpoints[(m / k) * hs..(m / k + 1) * hs].copy_from_slice(&cur);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    drop((prev, cur));

    // Backward, one block at a time.
    let group = panel.group_all().to_vec();
    let sizes = panel.sizes().to_vec();
    let mut buf = vec![0.0; k * hs];
    let mut loaded = usize::MAX;
    let mut beta = vec![1.0 / hs as f64; hs];
    let mut posterior = vec![0.0; sites * groups];
    let mut post = vec![0.0; groups];
    let mut s_alpha = vec![0.0; groups];
    let mut g_beta = vec![0.0; groups];
    // The interval to the right of the current site, once its backward step
    // is done: its K, its length and the scale beta was divided by. Its jump
    // terms and its G stay in `step` and `g_beta` until the next backward
    // step overwrites them, after they are read.
    let mut pending: Option<(f64, f64, f64)> = None;
    let mut step = Step::new(groups);
    for m in (0..sites).rev() {
        let b = m / k;
        if b != loaded {
            let start = b * k;
            let len = k.min(sites - start);
            buf[..hs].copy_from_slice(&checkpoints[b * hs..(b + 1) * hs]);
            for i in 1..len {
                let (done, rest) = buf.split_at_mut(i * hs);
                e.forward(start + i, Some(&done[(i - 1) * hs..]), &mut rest[..hs]);
            }
            loaded = b;
        }
        let alpha = &buf[(m - b * k) * hs..(m - b * k + 1) * hs];

        // The posterior at m.
        post.iter_mut().for_each(|p| *p = 0.0);
        s_alpha.iter_mut().for_each(|s| *s = 0.0);
        let mut mism = 0.0;
        for h in 0..hs {
            let w = alpha[h] * beta[h];
            post[group[h]] += w;
            s_alpha[group[h]] += alpha[h];
            if panel.allele(m, h) != target.get(m) {
                mism += w;
            }
        }
        let tot: f64 = post.iter().sum();
        for (a, p) in post.iter_mut().enumerate() {
            *p /= tot;
            stats.ancestry[a] += *p;
        }
        posterior[m * groups..(m + 1) * groups].copy_from_slice(&post);
        stats.sites += 1.0;
        stats.mismatches += mism / tot;

        // The interval from m to m+1: its expected events and jumps.
        if let Some((kk, d, scale)) = pending.take() {
            // The sum over both ends of alpha_m P e beta_{m+1}.
            let z = scale * tot;
            stats.events += kk / z;
            for a in 0..groups {
                stats.jumps[a] += step.jump[a] * s_alpha[a] * g_beta[a] / z;
                stats.exposure[a] += d * post[a];
            }
            stats.morgans += d;
        }

        // The backward step from m to m-1.
        if m > 0 {
            let d = (cm[m] - cm[m - 1]) / 100.0;
            step.set(params, &sizes, d);
            g_beta.iter_mut().for_each(|g| *g = 0.0);
            for (h, bh) in beta.iter_mut().enumerate() {
                *bh *= e.emission(m, h);
                g_beta[group[h]] += *bh;
            }
            let kk: f64 = (0..groups).map(|a| step.enter[a] * g_beta[a]).sum();
            for (h, bh) in beta.iter_mut().enumerate() {
                let a = group[h];
                *bh = step.stay[a] * *bh + step.jump[a] * g_beta[a] + kk;
            }
            let scale: f64 = beta.iter().sum();
            beta.iter_mut().for_each(|bh| *bh /= scale);
            pending = Some((kk, d, scale));
        }
    }
    Fb { posterior, stats }
}
