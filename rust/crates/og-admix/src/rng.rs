//! A seeded random number generator of our own: xoshiro256** (Blackman and
//! Vigna), its state filled by SplitMix64 from a 64-bit seed, as the
//! authors' reference code seeds it. No dependency, and the same seed gives
//! the same numbers on every machine and every build, which is what a
//! simulation with known truth and a bootstrap with fixed seeds need.

/// xoshiro256**: 256 bits of state, period 2^256 - 1.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rng {
    s: [u64; 4],
}

/// SplitMix64's step: the next output, advancing `state`.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// FNV-1a, 64 bits: a stable hash of a label, the same on every build.
fn fnv1a(label: &str) -> u64 {
    label.bytes().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

impl Rng {
    /// The generator for `seed`: four SplitMix64 outputs from it.
    pub fn new(seed: u64) -> Rng {
        let mut sm = seed;
        Rng {
            s: [
                splitmix64(&mut sm),
                splitmix64(&mut sm),
                splitmix64(&mut sm),
                splitmix64(&mut sm),
            ],
        }
    }

    /// An independent stream for one part of a job (a chromosome and
    /// haplotype, a bootstrap replicate), so a part's numbers do not change
    /// when other parts are added, removed or run in another order.
    pub fn stream(seed: u64, label: &str) -> Rng {
        Rng::new(seed ^ fnv1a(label))
    }

    /// The next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// A uniform number in [0, 1), from the top 53 bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// A uniform whole number in [0, n), without bias (Lemire's method).
    /// `n` must be at least 1.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "below(0) has no value to give");
        let threshold = n.wrapping_neg() % n;
        loop {
            let m = u128::from(self.next_u64()) * u128::from(n);
            if (m as u64) >= threshold {
                return (m >> 64) as u64;
            }
        }
    }

    /// An exponential waiting time with this `rate` (mean `1 / rate`);
    /// infinite for a rate of 0 or less, so an event that never happens
    /// never comes.
    pub fn exponential(&mut self, rate: f64) -> f64 {
        if rate <= 0.0 {
            return f64::INFINITY;
        }
        // 1 - u is in (0, 1], so the log is finite.
        -(1.0 - self.uniform()).ln() / rate
    }

    /// An index drawn with probability proportional to its weight, from the
    /// running totals of the weights (`cumulative[i]` is the sum of weights
    /// 0..=i; the last one is the total, above 0).
    pub fn pick(&mut self, cumulative: &[f64]) -> usize {
        let total = *cumulative.last().expect("at least one weight");
        let u = self.uniform() * total;
        cumulative
            .partition_point(|&c| c <= u)
            .min(cumulative.len() - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first five outputs of the authors' reference C code
    /// (xoshiro256starstar.c, state from splitmix64.c), compiled and run
    /// for these seeds.
    #[test]
    fn matches_the_reference_code() {
        let seed0 = if crate::PLANTED { 1 } else { 0 };
        let mut r = Rng::new(seed0);
        let got: Vec<u64> = (0..5).map(|_| r.next_u64()).collect();
        assert_eq!(
            got,
            [
                0x99ec_5f36_cb75_f2b4,
                0xbf6e_1f78_4956_452a,
                0x1a5f_849d_4933_e6e0,
                0x6aa5_94f1_262d_2d2c,
                0xbba5_ad4a_1f84_2e59
            ]
        );
        let mut r = Rng::new(42);
        let got: Vec<u64> = (0..5).map(|_| r.next_u64()).collect();
        assert_eq!(
            got,
            [
                0x1578_0b2e_0c2e_c716,
                0x6104_d986_6d11_3a7e,
                0xae17_5332_39e4_99a1,
                0xecb8_ad47_03b3_60a1,
                0xfde6_dc7f_e2ec_5e64
            ]
        );
        // The bare step, from the state 1, 2, 3, 4.
        let mut r = Rng { s: [1, 2, 3, 4] };
        let got: Vec<u64> = (0..5).map(|_| r.next_u64()).collect();
        assert_eq!(
            got,
            [
                0x2d00,
                0,
                0x5a00_7080,
                0x10e0_0000_0000_9d80,
                0x10e0_b61c_e100_9d80
            ]
        );
    }

    #[test]
    fn streams_are_stable_and_distinct() {
        let a: Vec<u64> = {
            let mut r = Rng::stream(7, "chr22/0");
            (0..4).map(|_| r.next_u64()).collect()
        };
        let b: Vec<u64> = {
            let mut r = Rng::stream(7, "chr22/0");
            (0..4).map(|_| r.next_u64()).collect()
        };
        let c: Vec<u64> = {
            let mut r = Rng::stream(7, "chr22/1");
            (0..4).map(|_| r.next_u64()).collect()
        };
        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_eq!(fnv1a(""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a("a"), 0xaf63_dc4c_8601_ec8c);
    }

    #[test]
    fn draws_have_their_laws() {
        let mut r = Rng::new(3);
        let n = 200_000;
        // uniform: in [0, 1), mean 1/2.
        let mut sum = 0.0;
        for _ in 0..n {
            let u = r.uniform();
            assert!((0.0..1.0).contains(&u));
            sum += u;
        }
        assert!((sum / f64::from(n) - 0.5).abs() < 0.005);
        // below: every value, none outside, about equally often.
        let mut counts = [0u32; 7];
        for _ in 0..n {
            counts[r.below(7) as usize] += 1;
        }
        let want = f64::from(n) / 7.0;
        for c in counts {
            assert!((f64::from(c) - want).abs() < 0.03 * want, "{counts:?}");
        }
        assert_eq!(r.below(1), 0);
        // exponential: mean 1 / rate, and P(X > mean) = 1/e.
        let rate = 4.0;
        let draws: Vec<f64> = (0..n).map(|_| r.exponential(rate)).collect();
        let mean = draws.iter().sum::<f64>() / f64::from(n);
        assert!((mean * rate - 1.0).abs() < 0.01, "{mean}");
        let above = draws.iter().filter(|&&x| x > 1.0 / rate).count() as f64 / f64::from(n);
        assert!((above - (-1.0f64).exp()).abs() < 0.005, "{above}");
        assert_eq!(r.exponential(0.0), f64::INFINITY);
        // pick: in proportion to the weights.
        let cum = [0.2, 0.2, 1.0];
        let mut got = [0u32; 3];
        for _ in 0..n {
            got[r.pick(&cum)] += 1;
        }
        assert_eq!(got[1], 0, "a zero weight is never picked");
        assert!((f64::from(got[0]) / f64::from(n) - 0.2).abs() < 0.005);
    }
}
