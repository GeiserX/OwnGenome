//! Ancestry fractions projected on fixed group frequencies, in place of
//! `admixture -P` (step 26b).
//!
//! **The model.** At SNV j a person carries g_j copies (0, 1 or 2) of the
//! allele the frequencies count, and group k carries it at frequency
//! f_jk. With ancestry fractions q on the simplex (q_k >= 0, sum 1), each
//! copy is the allele with probability p_j = sum_k q_k f_jk, so the
//! log-likelihood is
//!
//! ```text
//! L(q) = sum_j w_j [ g_j ln p_j + (2 - g_j) ln(1 - p_j) ]
//! ```
//!
//! (w_j is how many times a bootstrap replicate drew SNV j, 1 otherwise;
//! the binomial constant is left out, as ADMIXTURE leaves it out). Each
//! p_j is linear in q and ln is concave, so L is concave on the simplex:
//! every local maximum is the global one, and with frequencies whose
//! columns are not linearly dependent the maximum is unique. Any solver
//! that converges finds it, which is why the answer can match ADMIXTURE's
//! projection up to where ADMIXTURE stops.
//!
//! **The solver.** The EM update for q with f fixed is
//! `q_k <- q_k r_k`, where
//!
//! ```text
//! r_k = sum_j w_j [ g_j f_jk / p_j + (2 - g_j)(1 - f_jk) / (1 - p_j) ] / (2 sum_j w_j)
//! ```
//!
//! is the gradient of L over 2 sum_j w_j. It keeps q on the simplex and
//! never lowers L. SQUAREM (Varadhan and Roland 2008, the SqS3 step
//! length) extrapolates from two EM steps; the step is shortened toward
//! the plain EM step until every fraction stays above 0, and the
//! extrapolated point is kept only when its L is at least that of the
//! first EM step, so L never goes down. The run stops when one EM step
//! moves no fraction by more than [`TOLERANCE`].
//!
//! **The check.** At the maximum the Karush-Kuhn-Tucker conditions hold:
//! r_k = 1 for every group with q_k > 0 and r_k <= 1 for every group with
//! q_k = 0. [`project`] measures both at the answer (the largest
//! q_k |r_k - 1| and the largest r_k - 1) and returns an error when either
//! is above its tolerance, so a run cut short, or one stuck short of the
//! maximum, is never reported as an answer.
//!
//! **The bootstrap.** [`blocks`] cuts the SNVs into blocks of
//! [`BLOCK_CM`] on the genetic map, per chromosome, so no block holds two
//! chromosomes. A replicate draws as many blocks as there are, with
//! replacement ([`draw_weights`]), and draws each group's reference people
//! again ([`draw_people`]) to count its frequencies again, so the range
//! also holds the uncertainty of the frequencies themselves.

use crate::rng::Rng;
use crate::{Error, Result};
use std::collections::BTreeMap;

/// The run stops when an EM step moves no fraction by more than this.
pub const TOLERANCE: f64 = 1e-9;
/// The most passes over the SNVs before the run gives up.
pub const MAX_ITERATIONS: usize = 100_000;
/// The largest q_k |r_k - 1| the answer may have (complementary slackness).
pub const SLACK_TOLERANCE: f64 = 1e-7;
/// The largest r_k - 1 the answer may have (no group would gain).
pub const GAIN_TOLERANCE: f64 = 1e-4;
/// The length of a bootstrap block on the genetic map, in cM.
pub const BLOCK_CM: f64 = 5.0;

/// One projection: K groups, and per SNV the copies of the counted allele,
/// the K frequencies (row-major, SNV by group) and an optional weight.
#[derive(Debug, Clone, Copy)]
pub struct Problem<'a> {
    pub k: usize,
    pub genotypes: &'a [u8],
    pub freqs: &'a [f64],
    pub weights: Option<&'a [u32]>,
}

/// How far a run may go.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            tolerance: TOLERANCE,
            max_iterations: MAX_ITERATIONS,
        }
    }
}

/// The answer of a projection.
#[derive(Debug, Clone, PartialEq)]
pub struct Fit {
    /// The fractions, one per group, summing to 1.
    pub q: Vec<f64>,
    /// L at `q`.
    pub loglik: f64,
    /// Passes over the SNVs.
    pub iterations: usize,
    /// The largest q_k |r_k - 1| at `q`.
    pub slack: f64,
    /// The largest r_k - 1 at `q`.
    pub gain: f64,
}

impl<'a> Problem<'a> {
    /// Checks the shapes and the values: genotypes 0 to 2, every frequency
    /// strictly between 0 and 1, at least one SNV with weight.
    pub fn new(
        k: usize,
        genotypes: &'a [u8],
        freqs: &'a [f64],
        weights: Option<&'a [u32]>,
    ) -> Result<Problem<'a>> {
        if k == 0 {
            return Err(Error("a projection needs at least one group".into()));
        }
        if freqs.len() != genotypes.len() * k {
            return Err(Error(format!(
                "{} frequencies for {} SNVs of {k} groups",
                freqs.len(),
                genotypes.len()
            )));
        }
        if let Some(w) = weights
            && w.len() != genotypes.len()
        {
            return Err(Error(format!(
                "{} weights for {} SNVs",
                w.len(),
                genotypes.len()
            )));
        }
        if let Some(j) = genotypes.iter().position(|&g| g > 2) {
            return Err(Error(format!(
                "SNV {}: {} copies; a genotype holds 0, 1 or 2",
                j + 1,
                genotypes[j]
            )));
        }
        if let Some(i) = freqs.iter().position(|&f| !(f > 0.0 && f < 1.0)) {
            return Err(Error(format!(
                "SNV {}: frequency {} of group {} is not strictly between 0 and 1",
                i / k + 1,
                freqs[i],
                i % k + 1
            )));
        }
        let p = Problem {
            k,
            genotypes,
            freqs,
            weights,
        };
        if p.total() == 0.0 {
            return Err(Error("a projection needs at least one SNV".into()));
        }
        Ok(p)
    }

    fn weight(&self, j: usize) -> f64 {
        self.weights.map_or(1.0, |w| f64::from(w[j]))
    }

    /// The copies counted, 2 sum_j w_j.
    fn total(&self) -> f64 {
        2.0 * (0..self.genotypes.len())
            .map(|j| self.weight(j))
            .sum::<f64>()
    }

    /// L at `q`, and r (the gradient over 2 sum w) in `ratio`.
    pub fn evaluate(&self, q: &[f64], ratio: &mut [f64]) -> f64 {
        let k = self.k;
        ratio.iter_mut().for_each(|r| *r = 0.0);
        let (mut base, mut ll, mut total) = (0.0, 0.0, 0.0);
        for (j, &g) in self.genotypes.iter().enumerate() {
            let w = self.weight(j);
            if w == 0.0 {
                continue;
            }
            let row = &self.freqs[j * k..j * k + k];
            let p: f64 = q.iter().zip(row).map(|(a, b)| a * b).sum();
            let g = f64::from(g);
            let (mut a, mut b) = (0.0, 0.0);
            if g > 0.0 {
                a = g / p;
                ll += w * g * p.ln();
            }
            if g < 2.0 {
                b = (2.0 - g) / (1.0 - p);
                ll += w * (2.0 - g) * (1.0 - p).ln();
            }
            let d = w * (a - b);
            for (r, f) in ratio.iter_mut().zip(row) {
                *r += d * f;
            }
            base += w * b;
            total += 2.0 * w;
        }
        for r in ratio.iter_mut() {
            *r = (*r + base) / total;
        }
        ll
    }

    /// One EM step from `q` into `next`; L at `q`.
    fn em(&self, q: &[f64], next: &mut [f64], ratio: &mut [f64]) -> f64 {
        let ll = self.evaluate(q, ratio);
        for ((n, a), r) in next.iter_mut().zip(q).zip(ratio.iter()) {
            *n = a * r;
        }
        normalise(next);
        ll
    }
}

/// Scales `q` to sum 1.
fn normalise(q: &mut [f64]) {
    let s: f64 = q.iter().sum();
    if s > 0.0 {
        q.iter_mut().for_each(|x| *x /= s);
    }
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// The KKT measures at `q`: (largest q_k |r_k - 1|, largest r_k - 1).
pub fn kkt(q: &[f64], ratio: &[f64]) -> (f64, f64) {
    let slack = q
        .iter()
        .zip(ratio)
        .map(|(a, r)| a * (r - 1.0).abs())
        .fold(0.0, f64::max);
    let gain = ratio
        .iter()
        .map(|r| r - 1.0)
        .fold(f64::NEG_INFINITY, f64::max);
    (slack, gain)
}

/// The fractions that maximise L, from equal fractions, by EM with SQUAREM;
/// an error when the fit is not finite or the run ends before the KKT
/// conditions hold.
pub fn project(p: &Problem<'_>, opts: &Options) -> Result<Fit> {
    let k = p.k;
    let mut q = vec![1.0 / k as f64; k];
    let (mut q1, mut q2, mut qx, mut q3) = (vec![0.0; k], vec![0.0; k], vec![0.0; k], vec![0.0; k]);
    let mut ratio = vec![0.0; k];
    let (mut r, mut v) = (vec![0.0; k], vec![0.0; k]);
    let mut evals = 0usize;
    let mut converged = false;
    while evals < opts.max_iterations {
        p.em(&q, &mut q1, &mut ratio);
        evals += 1;
        if max_diff(&q1, &q) < opts.tolerance {
            q.copy_from_slice(&q1);
            converged = true;
            break;
        }
        if evals >= opts.max_iterations {
            q.copy_from_slice(&q1);
            break;
        }
        let l1 = p.em(&q1, &mut q2, &mut ratio);
        evals += 1;
        for i in 0..k {
            r[i] = q1[i] - q[i];
            v[i] = q2[i] - q1[i] - r[i];
        }
        let (rr, vv) = (dot(&r, &r), dot(&v, &v));
        let mut alpha = if vv > 0.0 { -(rr / vv).sqrt() } else { -1.0 };
        alpha = alpha.min(-1.0);
        // Shorten the step toward the plain EM step (alpha -1, which is
        // q2) until every fraction stays above 0.
        loop {
            if alpha > -1.0 - 1e-6 {
                qx.copy_from_slice(&q2);
                break;
            }
            for i in 0..k {
                qx[i] = q[i] - 2.0 * alpha * r[i] + alpha * alpha * v[i];
            }
            if qx.iter().all(|&x| x > 0.0) {
                normalise(&mut qx);
                break;
            }
            alpha = (alpha - 1.0) / 2.0;
        }
        if evals >= opts.max_iterations {
            q.copy_from_slice(&q2);
            break;
        }
        let lx = p.em(&qx, &mut q3, &mut ratio);
        evals += 1;
        // EM never lowers L, so EM(qx) is at least L(qx); q2 is at least
        // L(q1). Either way L does not go down.
        if lx.is_finite() && lx >= l1 {
            q.copy_from_slice(&q3);
        } else {
            q.copy_from_slice(&q2);
        }
    }
    let loglik = p.evaluate(&q, &mut ratio);
    let (slack, gain) = kkt(&q, &ratio);
    let fit = Fit {
        q,
        loglik,
        iterations: evals,
        slack,
        gain,
    };
    // A NaN passes every comparison below (and f64::max drops it), so a
    // fit that is not finite is refused first.
    if !(fit.q.iter().all(|x| x.is_finite())
        && fit.loglik.is_finite()
        && slack.is_finite()
        && gain.is_finite())
    {
        return Err(Error(format!(
            "the projection is not a finite number after {} passes (log-likelihood {}, fractions {:?}): a group frequency is too close to 0 or 1",
            fit.iterations, fit.loglik, fit.q
        )));
    }
    if !converged || slack > SLACK_TOLERANCE || gain > GAIN_TOLERANCE {
        return Err(Error(format!(
            "the projection stopped after {} passes short of the maximum: the KKT conditions do not hold (largest q|r - 1| {:e}, largest r - 1 {:e}{})",
            fit.iterations,
            fit.slack,
            fit.gain,
            if converged {
                String::new()
            } else {
                format!(
                    "; no step below {:e} within {} passes",
                    opts.tolerance, opts.max_iterations
                )
            }
        )));
    }
    Ok(fit)
}

/// The bootstrap blocks: the SNVs grouped by chromosome and by
/// `floor(cM / size)`, so no block holds two chromosomes; `chrom[j]` is any
/// number that names SNV j's chromosome. Blocks come in (chromosome, cM)
/// order, each with its SNVs in input order; empty blocks do not exist.
pub fn blocks(chrom: &[u32], cm: &[f64], size: f64) -> Vec<Vec<usize>> {
    let mut by: BTreeMap<(u32, i64), Vec<usize>> = BTreeMap::new();
    for (j, (&c, &x)) in chrom.iter().zip(cm).enumerate() {
        by.entry((c, (x / size).floor() as i64))
            .or_default()
            .push(j);
    }
    by.into_values().collect()
}

/// A replicate's weight of each of `snvs` SNVs: as many block draws, with
/// replacement, as there are blocks.
pub fn draw_weights(blocks: &[Vec<usize>], snvs: usize, rng: &mut Rng) -> Vec<u32> {
    let mut w = vec![0u32; snvs];
    if blocks.is_empty() {
        return w;
    }
    for _ in 0..blocks.len() {
        let b = rng.below(blocks.len() as u64) as usize;
        for &j in &blocks[b] {
            w[j] += 1;
        }
    }
    w
}

/// Each of `people` people's draw count when each group's members (person
/// indices) are drawn again, with replacement, as many times as the group
/// has members. A person in no group is drawn 0 times.
pub fn draw_people(groups: &[Vec<usize>], people: usize, rng: &mut Rng) -> Vec<u32> {
    let mut w = vec![0u32; people];
    for members in groups {
        for _ in 0..members.len() {
            w[members[rng.below(members.len() as u64) as usize]] += 1;
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A simulated K-group problem with fractions `q` and `snvs` SNVs.
    /// Each group's frequency is 0.02 or 0.98, drawn at random, and the
    /// copies follow p by error diffusion per frequency pattern instead of
    /// independent draws: each pattern's allele count is then 2 n p rounded,
    /// so the planted q is the maximum up to rounding, and the test measures
    /// the solver, not the luck of the draw.
    fn simulate(q: &[f64], snvs: usize, seed: u64) -> (Vec<u8>, Vec<f64>) {
        let k = q.len();
        let mut rng = Rng::new(seed);
        let mut acc = vec![0.0f64; 1 << k];
        let mut f = Vec::with_capacity(snvs * k);
        let mut g = Vec::with_capacity(snvs);
        for _ in 0..snvs {
            let bits = rng.below(1 << k) as usize;
            let row: Vec<f64> = (0..k)
                .map(|i| if bits >> i & 1 == 1 { 0.98 } else { 0.02 })
                .collect();
            let p: f64 = q.iter().zip(&row).map(|(a, b)| a * b).sum();
            let before = acc[bits].floor();
            acc[bits] += 2.0 * p;
            g.push((acc[bits].floor() - before) as u8);
            f.extend(row);
        }
        (g, f)
    }

    #[test]
    fn recovers_planted_fractions() {
        // Red-first truth test: three groups, 5,000 SNVs, the planted q
        // within 0.01. The red control swaps two frequency columns, so the
        // answer comes out with those two fractions swapped and fails.
        let truth = [0.6, 0.3, 0.1];
        let (g, mut f) = simulate(&truth, 5_000, 7);
        if crate::PLANTED {
            for row in f.chunks_mut(3) {
                row.swap(0, 1);
            }
        }
        let p = Problem::new(3, &g, &f, None).unwrap();
        let fit = project(&p, &Options::default()).unwrap();
        for (got, want) in fit.q.iter().zip(truth) {
            assert!((got - want).abs() <= 0.01, "q {:?} vs {truth:?}", fit.q);
        }
        assert!((fit.q.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(fit.slack <= SLACK_TOLERANCE && fit.gain <= GAIN_TOLERANCE);
    }

    #[test]
    fn a_truncated_run_fails_the_kkt_check() {
        let (g, f) = simulate(&[0.6, 0.3, 0.1], 5_000, 7);
        let p = Problem::new(3, &g, &f, None).unwrap();
        let opts = Options {
            max_iterations: 3,
            ..Options::default()
        };
        let e = project(&p, &opts).unwrap_err().0;
        assert!(
            e.contains("KKT conditions do not hold") && e.contains("within 3 passes"),
            "{e}"
        );
        // A plain EM run of the same length is far from the maximum too.
        let mut q = vec![1.0 / 3.0; 3];
        let (mut next, mut ratio) = (vec![0.0; 3], vec![0.0; 3]);
        for _ in 0..3 {
            p.em(&q.clone(), &mut next, &mut ratio);
            q.copy_from_slice(&next);
        }
        p.evaluate(&q, &mut ratio);
        assert!(kkt(&q, &ratio).0 > SLACK_TOLERANCE);
    }

    #[test]
    fn a_loose_stop_fails_the_kkt_check() {
        // A run that stops by its step size, well short of the maximum
        // pass count, can still be short of the optimum: only the KKT
        // measures catch it, not the `converged` flag.
        let (g, f) = simulate(&[0.6, 0.3, 0.1], 5_000, 7);
        let p = Problem::new(3, &g, &f, None).unwrap();
        let opts = Options {
            tolerance: 1e-3,
            ..Options::default()
        };
        let e = project(&p, &opts).unwrap_err().0;
        assert!(
            e.contains("KKT conditions do not hold") && !e.contains("within"),
            "{e}"
        );
    }

    #[test]
    fn the_answer_is_the_maximum() {
        // L at the answer is at least L at every point of a grid on the
        // simplex, and a group absent from the person ends at 0 with no
        // gain from growing it.
        let truth = [0.75, 0.25, 0.0];
        let (g, f) = simulate(&truth, 4_000, 3);
        let p = Problem::new(3, &g, &f, None).unwrap();
        let fit = project(&p, &Options::default()).unwrap();
        let mut ratio = vec![0.0; 3];
        for a in 0..=20 {
            for b in 0..=(20 - a) {
                let q = [a as f64 / 20.0, b as f64 / 20.0, (20 - a - b) as f64 / 20.0];
                let l = p.evaluate(&q, &mut ratio);
                assert!(
                    fit.loglik >= l - 1e-9 * l.abs(),
                    "{q:?}: {l} > {}",
                    fit.loglik
                );
            }
        }
        assert!(fit.q[2] < 0.01, "{:?}", fit.q);
    }

    #[test]
    fn plain_em_and_squarem_agree() {
        let (g, f) = simulate(&[0.2, 0.5, 0.3], 3_000, 11);
        let p = Problem::new(3, &g, &f, None).unwrap();
        let fit = project(&p, &Options::default()).unwrap();
        let mut q = vec![1.0 / 3.0; 3];
        let (mut next, mut ratio) = (vec![0.0; 3], vec![0.0; 3]);
        let mut last = f64::NEG_INFINITY;
        for _ in 0..200_000 {
            let l = p.em(&q.clone(), &mut next, &mut ratio);
            assert!(l >= last - 1e-9 * l.abs(), "EM lowered L");
            last = l;
            if max_diff(&next, &q) < 1e-12 {
                break;
            }
            q.copy_from_slice(&next);
        }
        assert!(max_diff(&q, &fit.q) < 1e-6, "{q:?} vs {:?}", fit.q);
        assert!(
            fit.iterations < 2_000,
            "SQUAREM took {} passes",
            fit.iterations
        );
    }

    #[test]
    fn weights_count_as_repeated_snvs() {
        let (g, f) = simulate(&[0.5, 0.5], 400, 5);
        let w: Vec<u32> = (0..400).map(|j| (j % 3) as u32).collect();
        let (mut g2, mut f2) = (Vec::new(), Vec::new());
        for j in 0..400 {
            for _ in 0..w[j] {
                g2.push(g[j]);
                f2.extend_from_slice(&f[j * 2..j * 2 + 2]);
            }
        }
        let a = project(
            &Problem::new(2, &g, &f, Some(&w)).unwrap(),
            &Options::default(),
        )
        .unwrap();
        let b = project(
            &Problem::new(2, &g2, &f2, None).unwrap(),
            &Options::default(),
        )
        .unwrap();
        assert!(max_diff(&a.q, &b.q) < 1e-7);
        assert!((a.loglik - b.loglik).abs() < 1e-6 * b.loglik.abs());
    }

    #[test]
    fn bad_input_is_an_error() {
        assert!(Problem::new(2, &[0, 3], &[0.5; 4], None).is_err());
        assert!(Problem::new(2, &[0, 1], &[0.5, 0.5, 1.0, 0.5], None).is_err());
        assert!(Problem::new(2, &[0, 1], &[0.5; 3], None).is_err());
        assert!(Problem::new(2, &[0, 1], &[0.5; 4], Some(&[0, 0])).is_err());
        assert!(Problem::new(0, &[], &[], None).is_err());
    }

    #[test]
    fn a_nonfinite_fit_is_refused() {
        // A frequency of 1e-310 is strictly between 0 and 1, but 2 / 1e-310
        // overflows: the EM step turns q into NaN, and f64::max drops the
        // NaN, so the step size reads 0 and the run looks converged.
        let p = Problem::new(1, &[2], &[1e-310], None).unwrap();
        let e = project(&p, &Options::default()).unwrap_err();
        assert!(e.0.contains("not a finite number"), "{}", e.0);
    }

    /// Blocks cut by SNV count, the bootstrap before this one: the red
    /// control of the chromosome test.
    fn blocks_by_count(n: usize, per: usize) -> Vec<Vec<usize>> {
        (0..n)
            .collect::<Vec<_>>()
            .chunks(per)
            .map(<[usize]>::to_vec)
            .collect()
    }

    #[test]
    fn no_block_spans_two_chromosomes() {
        // Two chromosomes of 7 and 5 SNVs, 0 to 12 cM: blocks of 5 cM.
        let chrom = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
        let cm = [0.0, 1.0, 4.9, 5.0, 6.0, 11.0, 12.0, 0.5, 2.0, 7.0, 8.0, 9.9];
        let b = if crate::PLANTED {
            blocks_by_count(12, 4)
        } else {
            blocks(&chrom, &cm, BLOCK_CM)
        };
        for block in &b {
            let c = chrom[block[0]];
            assert!(
                block.iter().all(|&j| chrom[j] == c),
                "block {block:?} spans two chromosomes"
            );
        }
        assert_eq!(
            b,
            vec![
                vec![0, 1, 2],
                vec![3, 4],
                vec![5, 6],
                vec![7, 8],
                vec![9, 10, 11]
            ]
        );
        let mut rng = Rng::new(1);
        let w = draw_weights(&b, 12, &mut rng);
        assert!(w.iter().any(|&x| x > 0));
        // A replicate draws 5 blocks: its SNV count is the sum of 5 block sizes.
        let sizes: Vec<usize> = b.iter().map(Vec::len).collect();
        let total: u32 = w.iter().sum();
        assert!(
            (sizes.iter().min().unwrap() * 5..=sizes.iter().max().unwrap() * 5)
                .contains(&(total as usize))
        );
    }

    #[test]
    fn people_draws_stay_in_their_group() {
        let groups = vec![vec![0, 2, 4], vec![1, 3]];
        let mut rng = Rng::new(9);
        for _ in 0..50 {
            let w = draw_people(&groups, 6, &mut rng);
            assert_eq!(w[0] + w[2] + w[4], 3);
            assert_eq!(w[1] + w[3], 2);
            assert_eq!(w[5], 0);
        }
    }
}
