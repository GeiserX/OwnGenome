use super::*;

/// The Python's smoothing, line for line: take the minimum of
/// `(agree, end - start, index)` over the whole list at every merge. O(R^2).
fn naive(sites: &[(i64, u8)]) -> Vec<Run> {
    let mut runs = runs_of(sites);
    while runs.len() > 1 {
        let Some(i) = (0..runs.len())
            .filter(|&k| weak(&runs[k]))
            .min_by_key(|&k| (runs[k].agree, runs[k].end - runs[k].start, k))
        else {
            break;
        };
        let lo = i.saturating_sub(1);
        let hi = (i + 1).min(runs.len() - 1);
        let keep = if lo != i {
            runs[lo].label
        } else {
            runs[hi].label
        };
        let mid: Vec<&Run> = runs[lo..=hi].iter().collect();
        let merged = merge(&runs[lo], &mid, &runs[hi], |k| lo + k == i, keep);
        runs.splice(lo..=hi, [merged]);
    }
    runs.retain(|r| r.agree >= MIN_SITES);
    runs
}

/// A small fixed generator, so a failure names its seed.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % n
    }
}

/// Sites with runs of every strength: switches are rare or frequent, the
/// spacing small or large, and positions repeat now and then (so starts tie).
fn sites(seed: u64) -> Vec<(i64, u8)> {
    let mut g = Lcg(seed);
    let n = 1 + g.next(600) as usize;
    let switch = 1 + g.next(40);
    let step = [1, 1_000, 25_000, 60_000][g.next(4) as usize];
    let (mut pos, mut label) = (1 + g.next(1_000) as i64, g.next(2) as u8);
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        if g.next(switch) == 0 {
            label ^= 1;
        }
        out.push((pos, label));
        if g.next(10) != 0 {
            pos += 1 + g.next(step) as i64;
        }
    }
    out
}

#[test]
fn heap_smooth_matches_the_python_order() {
    for seed in 0..3_000 {
        let s = sites(seed);
        assert_eq!(smooth(&s), naive(&s), "seed {seed}");
    }
}

/// Two weak runs of equal strength: the first one in the block goes first.
/// Taking the later one leaves a different map.
#[test]
fn ties_go_to_the_first_run() {
    let mut s = Vec::new();
    let mut pos = 0;
    let mut push = |label: u8, n: usize, s: &mut Vec<(i64, u8)>| {
        for _ in 0..n {
            pos += 30_000;
            s.push((pos, label));
        }
    };
    push(0, 25, &mut s);
    push(1, 3, &mut s);
    push(0, 3, &mut s);
    push(1, 25, &mut s);
    let got = smooth(&s);
    assert_eq!(got, naive(&s));
    assert_eq!(
        got.iter()
            .map(|r| (r.label, r.agree, r.discordant))
            .collect::<Vec<_>>(),
        [(0, 28, 3), (1, 25, 0)]
    );
}

#[test]
fn nothing_weak_stays_as_is() {
    let s: Vec<(i64, u8)> = (0..40).map(|i| (i * 100_000, u8::from(i >= 20))).collect();
    let got = smooth(&s);
    assert_eq!(got.len(), 2);
    assert!(smooth(&[]).is_empty());
    assert!(smooth(&[(5, 1)]).is_empty());
}

#[test]
fn blocks_split_at_gaps() {
    let b = blocks([1, 2, MAX_GAP + 3, MAX_GAP + 4].into_iter());
    assert_eq!(b, [0..2, 2..4]);
}
