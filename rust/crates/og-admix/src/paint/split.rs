//! A decoupled two-level copying model.
//!
//! Level 1, per group a: a Li and Stephens copying HMM over the n_a
//! reference haplotypes of a alone (jump rate rho_a, mismatch epsilon),
//! run left to right and right to left. Its per-site predictive
//! probabilities P(x_m | x_<m, a) and P(x_m | x_>m, a) give the local
//! evidence ell_a(m), the mean of their logs. Each group's copying keeps
//! running whatever the ancestry is, so a within-group haplotype switch
//! never costs an ancestry switch.
//!
//! Level 2: an ancestry HMM over the A groups, events at rate T per Morgan
//! drawing the new ancestry from mu, whose emission at site m is
//! exp(lambda * ell_a(m)). lambda <= 1 tempers the evidence: the per-site
//! predictive terms of one copied segment are strongly correlated, so their
//! raw sum overstates the evidence for short stretches.

use super::Params;
use super::bits::{Bits, RefPanel};
use super::hmm::{Fb, Stats};

/// lambda, the tempering of the local evidence. Chosen on the ancestry
/// engines bench (tests/bench/ancestry-engines.sh): 0.15 and 0.2 pass it,
/// 0.3 paints a pure CSA genome as 10% WEA.
pub const LAMBDA: f64 = 0.2;

/// Site by site, the local evidence of each group: `ell[m * A + a]`.
pub fn local_evidence(panel: &RefPanel, cm: &[f64], target: &Bits, params: &Params) -> Vec<f64> {
    let sites = panel.sites();
    let na = panel.groups();
    let group = panel.group_all();
    let eps = params.epsilon;
    let mut ell = vec![0.0f64; sites * na];
    let x: Vec<u8> = (0..sites).map(|m| target.get(m)).collect();
    for a in 0..na {
        let idx: Vec<usize> = (0..panel.haplotypes()).filter(|&h| group[h] == a).collect();
        let n = idx.len() as f64;
        let rho = params.rho[a];
        let mut alpha = vec![0.0f64; idx.len()];
        // forward: left to right
        for m in 0..sites {
            let pr = if m == 0 {
                1.0
            } else {
                -(-rho * (cm[m] - cm[m - 1]) / 100.0).exp_m1()
            };
            let row = panel.row(m);
            let mut c = 0.0;
            for (i, &h) in idx.iter().enumerate() {
                let e = if ((row[h / 64] >> (h % 64)) & 1) as u8 == x[m] {
                    1.0 - eps
                } else {
                    eps
                };
                let v = e * ((1.0 - pr) * alpha[i] + pr / n);
                alpha[i] = v;
                c += v;
            }
            let inv = 1.0 / c;
            alpha.iter_mut().for_each(|v| *v *= inv);
            ell[m * na + a] = 0.5 * c.ln();
        }
        // backward: right to left, the same model run the other way
        for m in (0..sites).rev() {
            let pr = if m + 1 == sites {
                1.0
            } else {
                -(-rho * (cm[m + 1] - cm[m]) / 100.0).exp_m1()
            };
            let row = panel.row(m);
            let mut c = 0.0;
            for (i, &h) in idx.iter().enumerate() {
                let e = if ((row[h / 64] >> (h % 64)) & 1) as u8 == x[m] {
                    1.0 - eps
                } else {
                    eps
                };
                let v = e * ((1.0 - pr) * alpha[i] + pr / n);
                alpha[i] = v;
                c += v;
            }
            let inv = 1.0 / c;
            alpha.iter_mut().for_each(|v| *v *= inv);
            ell[m * na + a] += 0.5 * c.ln();
        }
    }
    ell
}

/// The ancestry HMM over `ell` (sites * A) with tempering `lambda`.
pub fn ancestry_fb(ell: &[f64], cm: &[f64], params: &Params, lambda: f64, na: usize) -> Fb {
    let sites = cm.len();
    let mu = &params.mu;
    let t = params.generations;
    let mut stats = Stats::zero(na);
    let mut em = vec![0.0f64; sites * na];
    let mut logoff = 0.0;
    for m in 0..sites {
        let r = &ell[m * na..(m + 1) * na];
        let mx = r.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        logoff += lambda * mx;
        for a in 0..na {
            em[m * na + a] = (lambda * (r[a] - mx)).exp();
        }
    }
    let mut alpha = vec![0.0f64; sites * na];
    for m in 0..sites {
        let mut c = 0.0;
        if m == 0 {
            for a in 0..na {
                let v = mu[a] * em[a];
                alpha[a] = v;
                c += v;
            }
        } else {
            let pt = -(-t * (cm[m] - cm[m - 1]) / 100.0).exp_m1();
            for a in 0..na {
                let v = ((1.0 - pt) * alpha[(m - 1) * na + a] + pt * mu[a]) * em[m * na + a];
                alpha[m * na + a] = v;
                c += v;
            }
        }
        for a in 0..na {
            alpha[m * na + a] /= c;
        }
        stats.loglik += c.ln();
    }
    stats.loglik += logoff;
    let mut posterior = vec![0.0f64; sites * na];
    let mut beta = vec![1.0f64; na];
    let mut g = vec![0.0f64; na];
    for m in (0..sites).rev() {
        let mut tot = 0.0;
        for a in 0..na {
            tot += alpha[m * na + a] * beta[a];
        }
        for a in 0..na {
            let p = alpha[m * na + a] * beta[a] / tot;
            posterior[m * na + a] = p;
            stats.ancestry[a] += p;
        }
        stats.sites += 1.0;
        if m == 0 {
            break;
        }
        let d = (cm[m] - cm[m - 1]) / 100.0;
        let pt = -(-t * d).exp_m1();
        for a in 0..na {
            g[a] = em[m * na + a] * beta[a];
        }
        let enter: f64 = (0..na).map(|a| mu[a] * g[a]).sum();
        let mut z = 0.0;
        let mut s = 0.0;
        for a in 0..na {
            let nb = (1.0 - pt) * g[a] + pt * enter;
            z += alpha[(m - 1) * na + a] * nb;
            s += alpha[(m - 1) * na + a];
            beta[a] = nb;
        }
        stats.events += s * pt * enter / z;
        stats.morgans += d;
        let sc: f64 = beta.iter().sum();
        beta.iter_mut().for_each(|b| *b /= sc);
    }
    Fb { posterior, stats }
}
