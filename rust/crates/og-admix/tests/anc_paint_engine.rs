//! Truth tests of the painting engine (`og_admix::paint`):
//!
//! - a two-group chromosome of 2,000 sites with one breakpoint is painted
//!   with the breakpoint within 0.5 cM and a bp accuracy of at least 0.98,
//!   and the same check fails when the emission is set to 0.5;
//! - forward-backward equals a brute-force O(H^2) model on 30 sites and 12
//!   haplotypes within 1e-9, and the checkpointed pass equals the one that
//!   keeps every site within 1e-12;
//! - the peak memory of a 22-chromosome painting stays under 1 GB, and the
//!   checkpoints keep a long chromosome far below what full storage takes.
//!   The 1 GB check is a smoke check: this genome peaks far under it, so it
//!   cannot go red at this size. The memory evidence is the long-chromosome
//!   test, whose planted-red run (full storage) does go red;
//! - two runs, on 1 and on 4 threads, give the same bytes.
//!
//! Memory is the resident memory Linux reports (/proc/self/status); on
//! another system those two checks only run the painting.
//!
//! The simulator of `og ancestry-fine simulate` writes unphased genotypes,
//! and the engine paints phased haplotypes, so these tests build their
//! phased haplotypes with the same RNG and map.

use og_admix::map::GeneticMap;
use og_admix::paint::{
    Bits, ChromInput, Options, Params, RefPanel, checkpoint_block, forward_backward, paint,
    write_params,
};
use og_admix::panel::words;
use og_admix::rng::Rng;
use std::sync::Mutex;

/// The tests run one at a time, so a peak belongs to one test.
static ONE_AT_A_TIME: Mutex<()> = Mutex::new(());

fn lock() -> std::sync::MutexGuard<'static, ()> {
    ONE_AT_A_TIME.lock().unwrap_or_else(|e| e.into_inner())
}

/// A line `<key>: N kB` of /proc/self/status, in bytes.
#[cfg(target_os = "linux")]
fn status_bytes(key: &str) -> usize {
    let text = std::fs::read_to_string("/proc/self/status").expect("/proc/self/status");
    let line = text
        .lines()
        .find(|l| l.starts_with(key))
        .unwrap_or_else(|| panic!("no {key} in /proc/self/status"));
    let kb: usize = line[key.len()..]
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()
        .expect("a number of kB");
    kb * 1024
}

/// The peak resident memory of the process while `f` runs, and how far it
/// rose above the resident memory before (Linux: VmHWM of
/// /proc/self/status, reset through /proc/self/clear_refs). `None` on a
/// system without them; CI runs on Linux, where they must work.
fn peak_of<T>(f: impl FnOnce() -> T) -> (T, Option<(usize, usize)>) {
    #[cfg(target_os = "linux")]
    {
        std::fs::write("/proc/self/clear_refs", "5").expect("reset the peak RSS");
        let before = status_bytes("VmRSS:");
        let out = f();
        let peak = status_bytes("VmHWM:");
        (out, Some((peak, peak.saturating_sub(before))))
    }
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!("no /proc on this system: memory not measured");
        (f(), None)
    }
}

const PLANTED: bool = cfg!(feature = "planted-red");

/// A haplotype made of pieces of `founders`, switching to another founder
/// at `switch_per_site` per site, with each allele flipped at `noise`.
fn mosaic(rng: &mut Rng, founders: &[Vec<u8>], switch_per_site: f64, noise: f64) -> Vec<u8> {
    let sites = founders[0].len();
    let mut f = rng.below(founders.len() as u64) as usize;
    (0..sites)
        .map(|m| {
            if rng.uniform() < switch_per_site {
                f = rng.below(founders.len() as u64) as usize;
            }
            founders[f][m] ^ u8::from(rng.uniform() < noise)
        })
        .collect()
}

/// A simulated population: `groups` groups of `founders` founder
/// haplotypes each, at `sites` sites, whose allele frequencies drift apart
/// from a shared one by up to 0.25.
fn founders(rng: &mut Rng, groups: usize, founders: usize, sites: usize) -> Vec<Vec<Vec<u8>>> {
    let base: Vec<f64> = (0..sites).map(|_| 0.1 + 0.8 * rng.uniform()).collect();
    (0..groups)
        .map(|_| {
            let freq: Vec<f64> = base
                .iter()
                .map(|p| (p + 0.5 * (rng.uniform() - 0.5)).clamp(0.02, 0.98))
                .collect();
            (0..founders)
                .map(|_| freq.iter().map(|&p| u8::from(rng.uniform() < p)).collect())
                .collect()
        })
        .collect()
}

/// A reference panel of `per_group` mosaic haplotypes of each group.
fn panel_of(rng: &mut Rng, pop: &[Vec<Vec<u8>>], per_group: usize) -> RefPanel {
    let sites = pop[0][0].len();
    let haps: Vec<(usize, Vec<u8>)> = pop
        .iter()
        .enumerate()
        .flat_map(|(g, f)| {
            (0..per_group)
                .map(|_| (g, mosaic(rng, f, 0.01, 0.005)))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut panel = RefPanel::new(haps.iter().map(|(g, _)| *g).collect(), pop.len()).unwrap();
    let mut row = vec![0u64; words(haps.len())];
    for m in 0..sites {
        row.iter_mut().for_each(|w| *w = 0);
        for (h, (_, a)) in haps.iter().enumerate() {
            row[h / 64] |= u64::from(a[m]) << (h % 64);
        }
        panel.push_site(&row);
    }
    panel
}

/// cM of each site of a chromosome of `sites` sites 1 kb apart, through
/// a map of 1 cM per Mb.
fn map_cm(sites: usize) -> Vec<f64> {
    let map = GeneticMap::from_markers(vec![(0, 0.0), (1_000_000_000, 1000.0)]).unwrap();
    (0..sites)
        .map(|m| map.cm_at(1_000 * (m as i64 + 1)))
        .collect()
}

// ---- Truth: one breakpoint ----

const SITES: usize = 2_000;
const BREAK: usize = 1_000;

struct Truth {
    chrom: ChromInput,
    /// The true group of each site, per haplotype.
    truth: [Vec<usize>; 2],
}

/// Two groups A (0) and B (1). Haplotype 0 copies A up to site 1,000 and B
/// after; haplotype 1 is B throughout.
fn one_breakpoint() -> Truth {
    let mut rng = Rng::stream(20261119, "anc-paint-engine/one-breakpoint");
    let pop = founders(&mut rng, 2, 4, SITES);
    let panel = panel_of(&mut rng, &pop, 30);
    let a = mosaic(&mut rng, &pop[0], 0.01, 0.005);
    let b = mosaic(&mut rng, &pop[1], 0.01, 0.005);
    let hap0: Vec<u8> = (0..SITES)
        .map(|m| if m < BREAK { a[m] } else { b[m] })
        .collect();
    let hap1 = mosaic(&mut rng, &pop[1], 0.01, 0.005);
    Truth {
        chrom: ChromInput {
            name: "chr1".into(),
            cm: map_cm(SITES),
            panel,
            targets: vec![Bits::from_alleles(&hap0), Bits::from_alleles(&hap1)],
        },
        truth: [
            (0..SITES).map(|m| usize::from(m >= BREAK)).collect(),
            vec![1; SITES],
        ],
    }
}

/// The fixture Q of the continental projection: a quarter A.
fn fixture_q() -> Vec<f64> {
    vec![0.25, 0.75]
}

/// Paints the truth chromosome from `start` with `rounds` EM rounds and
/// checks the result: bp accuracy at least 0.98 and a called breakpoint
/// within 0.5 cM of the true one.
fn paint_and_score(t: &Truth, start: Params, rounds: usize) -> Result<(), String> {
    let p = paint(
        std::slice::from_ref(&t.chrom),
        start,
        &Options { threads: 2, rounds },
    )
    .map_err(|e| e.0)?;
    let pos = |m: usize| 1_000 * (m as u64 + 1);
    let (mut right, mut all) = (0u64, 0u64);
    for (h, hap) in p.chroms[0].haplotypes.iter().enumerate() {
        for m in 0..SITES - 1 {
            let bp = pos(m + 1) - pos(m);
            all += bp;
            if hap.calls[m] == t.truth[h][m] {
                right += bp;
            }
        }
    }
    let accuracy = right as f64 / all as f64;
    let calls = &p.chroms[0].haplotypes[0].calls;
    let breaks: Vec<f64> = (1..SITES)
        .filter(|&m| calls[m] != calls[m - 1])
        .map(|m| (t.chrom.cm[m - 1] + t.chrom.cm[m]) / 2.0)
        .collect();
    let true_cm = (t.chrom.cm[BREAK - 1] + t.chrom.cm[BREAK]) / 2.0;
    let nearest = breaks
        .iter()
        .map(|b| (b - true_cm).abs())
        .fold(f64::INFINITY, f64::min);
    eprintln!(
        "bp accuracy {accuracy}, nearest called breakpoint {nearest} cM off, {} EM rounds, T {}, epsilon {}, mu {:?}",
        p.rounds, p.params.generations, p.params.epsilon, p.params.mu
    );
    if accuracy < 0.98 {
        return Err(format!("bp accuracy {accuracy}, below 0.98"));
    }
    if nearest > 0.5 {
        return Err(format!(
            "no called breakpoint within 0.5 cM of {true_cm} cM (called: {breaks:?})"
        ));
    }
    Ok(())
}

#[test]
fn paints_one_breakpoint_within_half_a_cm() {
    let _g = lock();
    let t = one_breakpoint();
    let mut start = Params::initial(
        fixture_q(),
        t.chrom.panel.sizes(),
        og_admix::paint::DEFAULT_NE,
    );
    if PLANTED {
        // The red control: an emission that says nothing.
        start.epsilon = 0.5;
    }
    let rounds = if PLANTED { 0 } else { 5 };
    paint_and_score(&t, start, rounds).unwrap();
}

#[test]
fn an_emission_of_one_half_fails_the_breakpoint_check() {
    let _g = lock();
    let t = one_breakpoint();
    let mut start = Params::initial(
        fixture_q(),
        t.chrom.panel.sizes(),
        og_admix::paint::DEFAULT_NE,
    );
    start.epsilon = 0.5;
    let got = paint_and_score(&t, start, 0);
    assert!(got.is_err(), "the check passed with no information");
}

// ---- Exactness: brute force and checkpoints ----

/// A small random case: groups of 5, 4 and 3 haplotypes, `sites` sites at
/// uneven distances (some 0), and parameters where every term matters.
fn small_case(seed: u64, sites: usize, sizes: &[usize]) -> (RefPanel, Vec<f64>, Bits, Params) {
    let mut rng = Rng::new(seed);
    let group: Vec<usize> = sizes
        .iter()
        .enumerate()
        .flat_map(|(g, &n)| std::iter::repeat_n(g, n))
        .collect();
    // Interleave the groups, so a group is not one run of haplotypes.
    let mut order: Vec<usize> = (0..group.len()).collect();
    for i in (1..order.len()).rev() {
        order.swap(i, rng.below(i as u64 + 1) as usize);
    }
    let group: Vec<usize> = order.iter().map(|&i| group[i]).collect();
    let hs = group.len();
    let mut panel = RefPanel::new(group, sizes.len()).unwrap();
    let mut row = vec![0u64; words(hs)];
    for _ in 0..sites {
        row.iter_mut().for_each(|w| *w = 0);
        for h in 0..hs {
            row[h / 64] |= u64::from(rng.uniform() < 0.4) << (h % 64);
        }
        panel.push_site(&row);
    }
    let mut cm = vec![0.0];
    for m in 1..sites {
        let step = if m % 7 == 3 { 0.0 } else { 1.2 * rng.uniform() };
        cm.push(cm[m - 1] + step);
    }
    let target = Bits::from_alleles(
        &(0..sites)
            .map(|_| u8::from(rng.uniform() < 0.5))
            .collect::<Vec<_>>(),
    );
    let mut mu: Vec<f64> = sizes.iter().map(|_| 0.2 + rng.uniform()).collect();
    let total: f64 = mu.iter().sum();
    mu.iter_mut().for_each(|x| *x /= total);
    let params = Params {
        mu,
        generations: 40.0,
        epsilon: 0.07,
        rho: sizes.iter().map(|_| 20.0 + 300.0 * rng.uniform()).collect(),
    };
    (panel, cm, target, params)
}

/// The same model with its full H x H transition matrix: the group
/// posteriors per site, the log-likelihood and the expected events and
/// jumps.
struct Brute {
    posterior: Vec<f64>,
    loglik: f64,
    events: f64,
    jumps: Vec<f64>,
}

// Written index by index, as the equations read.
#[allow(clippy::needless_range_loop)]
fn brute_force(panel: &RefPanel, cm: &[f64], target: &Bits, p: &Params) -> Brute {
    let hs = panel.haplotypes();
    let a_n = panel.groups();
    let sites = panel.sites();
    let g = |h: usize| panel.group_of(h);
    let n = |a: usize| panel.sizes()[a] as f64;
    let e = |m: usize, h: usize| {
        if panel.allele(m, h) == target.get(m) {
            1.0 - p.epsilon
        } else {
            p.epsilon
        }
    };
    // trans[m][i][j], with the event and jump parts apart.
    let parts = |m: usize| {
        let d = (cm[m] - cm[m - 1]) / 100.0;
        let pt = 1.0 - (-p.generations * d).exp();
        let mut full = vec![vec![0.0; hs]; hs];
        let mut event = vec![vec![0.0; hs]; hs];
        let mut jump = vec![vec![0.0; hs]; hs];
        for i in 0..hs {
            for j in 0..hs {
                let (a, b) = (g(i), g(j));
                let pr = 1.0 - (-p.rho[a] * d).exp();
                let ev = pt * p.mu[b] / n(b);
                let ju = if a == b { (1.0 - pt) * pr / n(a) } else { 0.0 };
                let st = if i == j { (1.0 - pt) * (1.0 - pr) } else { 0.0 };
                event[i][j] = ev;
                jump[i][j] = ju;
                full[i][j] = ev + ju + st;
            }
        }
        (full, event, jump)
    };
    let mut alpha = vec![vec![0.0; hs]; sites];
    let mut loglik = 0.0;
    for h in 0..hs {
        alpha[0][h] = p.mu[g(h)] / n(g(h)) * e(0, h);
    }
    let c: f64 = alpha[0].iter().sum();
    alpha[0].iter_mut().for_each(|x| *x /= c);
    loglik += c.ln();
    for m in 1..sites {
        let (full, _, _) = parts(m);
        for j in 0..hs {
            alpha[m][j] = (0..hs).map(|i| alpha[m - 1][i] * full[i][j]).sum::<f64>() * e(m, j);
        }
        let c: f64 = alpha[m].iter().sum();
        alpha[m].iter_mut().for_each(|x| *x /= c);
        loglik += c.ln();
    }
    let mut beta = vec![vec![1.0; hs]; sites];
    for m in (0..sites - 1).rev() {
        let (full, _, _) = parts(m + 1);
        for i in 0..hs {
            beta[m][i] = (0..hs)
                .map(|j| full[i][j] * e(m + 1, j) * beta[m + 1][j])
                .sum();
        }
        let c: f64 = beta[m].iter().sum();
        beta[m].iter_mut().for_each(|x| *x /= c);
    }
    let mut posterior = vec![0.0; sites * a_n];
    for m in 0..sites {
        let tot: f64 = (0..hs).map(|h| alpha[m][h] * beta[m][h]).sum();
        for h in 0..hs {
            posterior[m * a_n + g(h)] += alpha[m][h] * beta[m][h] / tot;
        }
    }
    let (mut events, mut jumps) = (0.0, vec![0.0; a_n]);
    for m in 1..sites {
        let (full, event, jump) = parts(m);
        let w =
            |i: usize, j: usize, t: &[Vec<f64>]| alpha[m - 1][i] * t[i][j] * e(m, j) * beta[m][j];
        let mut z = 0.0;
        for i in 0..hs {
            for j in 0..hs {
                z += w(i, j, &full);
            }
        }
        for i in 0..hs {
            for j in 0..hs {
                events += w(i, j, &event) / z;
                jumps[g(i)] += w(i, j, &jump) / z;
            }
        }
    }
    Brute {
        posterior,
        loglik,
        events,
        jumps,
    }
}

#[test]
fn forward_backward_equals_the_brute_force_model() {
    let _g = lock();
    for seed in [1u64, 2, 3] {
        let (panel, cm, target, params) = small_case(seed, 30, &[5, 4, 3]);
        assert_eq!(panel.haplotypes(), 12);
        let want = brute_force(&panel, &cm, &target, &params);
        let mut block = checkpoint_block(30);
        if PLANTED {
            // The red control: the engine with another model.
            block = 30;
        }
        let mut model = params.clone();
        if PLANTED {
            model.rho[0] *= 1.01;
        }
        let got = forward_backward(&panel, &cm, &target, &model, block);
        for (i, (x, y)) in got.posterior.iter().zip(&want.posterior).enumerate() {
            assert!(
                (x - y).abs() <= 1e-9,
                "seed {seed}: posterior {i}: {x} against the brute force {y}"
            );
        }
        assert!(
            (got.stats.loglik - want.loglik).abs() <= 1e-9,
            "seed {seed}: loglik {} against {}",
            got.stats.loglik,
            want.loglik
        );
        assert!(
            (got.stats.events - want.events).abs() <= 1e-9,
            "seed {seed}: events {} against {}",
            got.stats.events,
            want.events
        );
        for a in 0..3 {
            assert!(
                (got.stats.jumps[a] - want.jumps[a]).abs() <= 1e-9,
                "seed {seed}: jumps of group {a}: {} against {}",
                got.stats.jumps[a],
                want.jumps[a]
            );
        }
    }
}

#[test]
fn checkpoints_give_what_full_storage_gives() {
    let _g = lock();
    let (panel, cm, target, params) = small_case(7, 500, &[20, 12, 8]);
    let full = forward_backward(&panel, &cm, &target, &params, 500);
    for block in [checkpoint_block(500), 1, 7, 64, 499] {
        let block = if PLANTED { block + 1 } else { block };
        let got = forward_backward(&panel, &cm, &target, &params, block);
        let mut model = params.clone();
        if PLANTED {
            model.epsilon += 1e-6;
        }
        let got = if PLANTED {
            forward_backward(&panel, &cm, &target, &model, block)
        } else {
            got
        };
        for (i, (x, y)) in got.posterior.iter().zip(&full.posterior).enumerate() {
            assert!(
                (x - y).abs() <= 1e-12,
                "block {block}: posterior {i}: {x} against {y}"
            );
        }
        assert!((got.stats.loglik - full.stats.loglik).abs() <= 1e-12);
        assert!((got.stats.events - full.stats.events).abs() <= 1e-12);
        assert!((got.stats.mismatches - full.stats.mismatches).abs() <= 1e-12);
    }
}

// ---- Memory and determinism ----

/// A genome of 22 chromosomes, three groups of 40 reference haplotypes and
/// one target person.
fn genome() -> Vec<ChromInput> {
    (1..=22)
        .map(|c| {
            let mut rng = Rng::stream(20261119, &format!("anc-paint-engine/chr{c}"));
            let sites = 600 + 20 * c;
            let pop = founders(&mut rng, 3, 4, sites);
            let panel = panel_of(&mut rng, &pop, 40);
            let a = mosaic(&mut rng, &pop[c % 3], 0.01, 0.005);
            let b = mosaic(&mut rng, &pop[(c + 1) % 3], 0.01, 0.005);
            let cut = sites / 3;
            let hap0: Vec<u8> = (0..sites)
                .map(|m| if m < cut { a[m] } else { b[m] })
                .collect();
            let hap1 = mosaic(&mut rng, &pop[(c + 2) % 3], 0.01, 0.005);
            ChromInput {
                name: format!("chr{c}"),
                cm: map_cm(sites),
                panel,
                targets: vec![Bits::from_alleles(&hap0), Bits::from_alleles(&hap1)],
            }
        })
        .collect()
}

/// Every byte of a painting: params.tsv, then each haplotype's posteriors
/// and calls.
fn bytes(chroms: &[ChromInput], threads: usize) -> Vec<u8> {
    let start = Params::initial(
        vec![0.3, 0.3, 0.4],
        chroms[0].panel.sizes(),
        og_admix::paint::DEFAULT_NE,
    );
    let p = paint(chroms, start, &Options { threads, rounds: 5 }).unwrap();
    let mut out = Vec::new();
    write_params(&mut out, &p, &["A".into(), "B".into(), "C".into()]).unwrap();
    for c in &p.chroms {
        out.extend_from_slice(c.name.as_bytes());
        for h in &c.haplotypes {
            for x in &h.posterior {
                out.extend_from_slice(&x.to_le_bytes());
            }
            out.extend(h.calls.iter().map(|&a| a as u8));
        }
    }
    out
}

/// The 1 GB bound is a smoke check (this genome peaks at a few MB); the
/// memory evidence is [`checkpoints_keep_a_long_chromosome_small`].
#[test]
fn a_22_chromosome_painting_stays_under_1_gb_and_repeats_byte_for_byte() {
    let _g = lock();
    let chroms = genome();
    let (one, peak) = peak_of(|| bytes(&chroms, 1));
    if let Some((peak, _)) = peak {
        eprintln!("peak resident memory {peak} bytes");
        assert!(peak < 1 << 30, "peak {peak} bytes, 1 GB or more");
    }
    let four = bytes(&chroms, 4);
    let mut again = bytes(&chroms, 4);
    if PLANTED {
        // The red control: one bit of one posterior differs.
        let last = again.len() - 1;
        again[last] ^= 1;
    }
    assert!(one == four, "1 and 4 threads gave different bytes");
    assert!(four == again, "two runs gave different bytes");
    let text =
        String::from_utf8_lossy(&one[..one.iter().position(|&b| b == b'c').unwrap()]).into_owned();
    assert!(
        text.starts_with("#parameter\tgroup\tvalue\nT\t.\t"),
        "{text}"
    );
}

#[test]
fn checkpoints_keep_a_long_chromosome_small() {
    let _g = lock();
    // 40,000 sites and 400 haplotypes: keeping every site's forward values
    // takes 40,000 x 400 x 8 = 128 MB.
    let (panel, cm, target, params) = small_case(11, 40_000, &[150, 150, 100]);
    let full_alpha = 40_000 * 400 * 8;
    let bound = full_alpha / 16;
    let block = if PLANTED {
        40_000
    } else {
        checkpoint_block(40_000)
    };
    let (_, peak) = peak_of(|| forward_backward(&panel, &cm, &target, &params, block));
    let Some((_, rise)) = peak else { return };
    assert!(
        rise < bound,
        "memory rose {rise} bytes with checkpoints, not below {bound}"
    );
    // The control: keeping every site goes over the bound.
    let (_, full) = peak_of(|| forward_backward(&panel, &cm, &target, &params, 40_000));
    let (_, full_rise) = full.expect("measured above");
    eprintln!("memory rose {rise} bytes with checkpoints, {full_rise} keeping every site");
    assert!(
        full_rise > bound,
        "full storage rose {full_rise} bytes, the bound {bound} would not see it"
    );
}

#[test]
fn bad_inputs_stop_with_a_reason() {
    let _g = lock();
    let t = one_breakpoint();
    let ok = Params::initial(fixture_q(), t.chrom.panel.sizes(), 10_000.0);
    let run = |c: &ChromInput, p: Params, threads: usize| {
        paint(std::slice::from_ref(c), p, &Options { threads, rounds: 0 })
            .map(|_| ())
            .map_err(|e| e.0)
    };
    let mut bad = ok.clone();
    bad.mu = vec![0.5, 0.6];
    assert!(run(&t.chrom, bad, 1).unwrap_err().contains("sum to"));
    let mut bad = ok.clone();
    bad.epsilon = 0.6;
    assert!(run(&t.chrom, bad, 1).unwrap_err().contains("epsilon"));
    let mut c = t.chrom.clone();
    c.cm[5] = c.cm[4] - 1.0;
    assert!(
        run(&c, ok.clone(), 1)
            .unwrap_err()
            .contains("never go down")
    );
    let mut c = t.chrom.clone();
    c.targets[1] = Bits::from_alleles(&[0, 1]);
    assert!(
        run(&c, ok.clone(), 1)
            .unwrap_err()
            .contains("target haplotype 1")
    );
    let mut c = t.chrom.clone();
    c.targets.clear();
    assert!(
        run(&c, ok.clone(), 1)
            .unwrap_err()
            .contains("no target haplotype")
    );
    assert!(run(&t.chrom, ok, 0).unwrap_err().contains("thread"));
    assert!(
        RefPanel::new(vec![0, 0, 2], 3)
            .unwrap_err()
            .0
            .contains("group 1")
    );
}
