//! Regression tests of the two-level painting engine (`og_admix::paint`)
//! against the faults of the joint (ancestry, haplotype) model it
//! replaced. Each test fails on the joint model at a3c86ee9:
//!
//! - an unadmixed haplotype from groups whose frequencies differ little is
//!   painted as one group: the joint model reads short matches to another
//!   group's haplotypes as ancestry switches, and EM then raises T;
//! - the fit keeps mu at the start (the continental projection); the joint
//!   model refits it from the posterior and drifts towards absent groups;
//! - a real minor ancestry is still painted: the two-level model does not
//!   buy its calm by erasing small groups.

use og_admix::map::GeneticMap;
use og_admix::paint::{Bits, ChromInput, Options, Params, RefPanel, paint};
use og_admix::panel::words;
use og_admix::rng::Rng;

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

/// `groups` groups of `founders` founder haplotypes at `sites` sites,
/// whose allele frequencies differ from a shared one by up to `drift`.
fn founders(
    rng: &mut Rng,
    groups: usize,
    founders: usize,
    sites: usize,
    drift: f64,
) -> Vec<Vec<Vec<u8>>> {
    let base: Vec<f64> = (0..sites).map(|_| 0.1 + 0.8 * rng.uniform()).collect();
    (0..groups)
        .map(|_| {
            let freq: Vec<f64> = base
                .iter()
                .map(|p| (p + 2.0 * drift * (rng.uniform() - 0.5)).clamp(0.02, 0.98))
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

/// cM of each site, 1 kb apart, through a map of 1 cM per Mb.
fn map_cm(sites: usize) -> Vec<f64> {
    let map = GeneticMap::from_markers(vec![(0, 0.0), (1_000_000_000, 1000.0)]).unwrap();
    (0..sites)
        .map(|m| map.cm_at(1_000 * (m as i64 + 1)))
        .collect()
}

const SITES: usize = 6_000;
const GROUPS: usize = 3;

/// Three close groups (frequencies within 0.1 of a shared one), 8
/// founders and 20 reference haplotypes each, and a target person whose
/// two haplotypes are fresh mosaics of group 0's founders, except
/// haplotype 1's sites in `minor`, copied from group 2.
fn world(seed: &str, minor: std::ops::Range<usize>) -> (ChromInput, [Vec<usize>; 2]) {
    let mut rng = Rng::stream(20261124, seed);
    let pop = founders(&mut rng, GROUPS, 8, SITES, 0.1);
    let panel = panel_of(&mut rng, &pop, 20);
    let hap0 = mosaic(&mut rng, &pop[0], 0.02, 0.005);
    let other = mosaic(&mut rng, &pop[2], 0.02, 0.005);
    let hap1: Vec<u8> = mosaic(&mut rng, &pop[0], 0.02, 0.005)
        .into_iter()
        .enumerate()
        .map(|(m, a)| if minor.contains(&m) { other[m] } else { a })
        .collect();
    let truth = [
        vec![0; SITES],
        (0..SITES)
            .map(|m| if minor.contains(&m) { 2 } else { 0 })
            .collect(),
    ];
    let chrom = ChromInput {
        name: "chr1".into(),
        cm: map_cm(SITES),
        panel,
        targets: vec![Bits::from_alleles(&hap0), Bits::from_alleles(&hap1)],
    };
    (chrom, truth)
}

/// Share of sites called right, and the number of called stretches, over
/// both haplotypes.
fn score(calls: [&[usize]; 2], truth: &[Vec<usize>; 2]) -> (f64, usize) {
    let mut right = 0;
    let mut stretches = 0;
    for (c, t) in calls.iter().zip(truth) {
        right += c.iter().zip(t).filter(|(a, b)| a == b).count();
        stretches += 1 + (1..c.len()).filter(|&m| c[m] != c[m - 1]).count();
    }
    (right as f64 / (2 * SITES) as f64, stretches)
}

fn start(chrom: &ChromInput, mu: Vec<f64>) -> Params {
    Params::initial(mu, chrom.panel.sizes(), og_admix::paint::DEFAULT_NE)
}

#[test]
fn an_unadmixed_person_is_painted_as_one_group() {
    let (chrom, truth) = world("unadmixed", 0..0);
    let mu = vec![0.98, 0.01, 0.01];
    let p = paint(
        std::slice::from_ref(&chrom),
        start(&chrom, mu),
        &Options {
            threads: 2,
            rounds: 5,
        },
    )
    .unwrap();
    let h = &p.chroms[0].haplotypes;
    let (accuracy, stretches) = score([&h[0].calls, &h[1].calls], &truth);
    eprintln!(
        "accuracy {accuracy}, {stretches} stretches (truth 2), T {}",
        p.params.generations
    );
    assert!(accuracy >= 0.99, "accuracy {accuracy}, below 0.99");
    assert!(stretches <= 4, "{stretches} stretches, truth 2");
    assert!(
        p.params.generations < 30.0,
        "T fitted to {}, from a start of 10 with no ancestry event",
        p.params.generations
    );
}

#[test]
fn the_fit_keeps_mu_at_the_projection() {
    let (chrom, _) = world("unadmixed", 0..0);
    let mu = vec![0.9, 0.05, 0.05];
    let p = paint(
        std::slice::from_ref(&chrom),
        start(&chrom, mu.clone()),
        &Options {
            threads: 2,
            rounds: 5,
        },
    )
    .unwrap();
    assert_eq!(p.params.mu, mu);
}

#[test]
fn a_minor_ancestry_of_four_percent_is_still_painted() {
    // Haplotype 1 carries group 2 on sites 2,000 to 2,480: 4% of the
    // person, 0.48 cM.
    let (chrom, truth) = world("minor", 2_000..2_480);
    let mu = vec![0.96, 0.0001, 0.0399];
    let p = paint(
        std::slice::from_ref(&chrom),
        start(&chrom, mu),
        &Options {
            threads: 2,
            rounds: 5,
        },
    )
    .unwrap();
    let h = &p.chroms[0].haplotypes;
    let minor = h[1].calls[2_000..2_480].iter().filter(|&&c| c == 2).count();
    let (accuracy, stretches) = score([&h[0].calls, &h[1].calls], &truth);
    eprintln!("minor sites called 2: {minor} of 480, accuracy {accuracy}, {stretches} stretches");
    assert!(
        minor >= 400,
        "{minor} of the 480 minor sites called as group 2"
    );
    assert!(accuracy >= 0.98, "accuracy {accuracy}, below 0.98");
}
