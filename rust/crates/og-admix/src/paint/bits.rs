//! The bit-packed inputs of the painting: the reference haplotypes of one
//! chromosome, one bit per haplotype per site, and a target haplotype, one
//! bit per site.

use crate::panel::{Haplotypes, bit, words};
use crate::{Error, Result};

/// The reference haplotypes of one chromosome and the group (ancestry) of
/// each. A site takes `ceil(H / 64)` words, so a panel of H haplotypes at M
/// sites holds `M * H / 8` bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefPanel {
    haplotypes: usize,
    words: usize,
    bits: Vec<u64>,
    group: Vec<usize>,
    size: Vec<usize>,
}

impl RefPanel {
    /// An empty panel: `group[h]` is the group of haplotype `h`, an index
    /// below `groups`. Every group needs at least one haplotype.
    pub fn new(group: Vec<usize>, groups: usize) -> Result<RefPanel> {
        if groups == 0 {
            return Err(Error("the painting needs at least one group".into()));
        }
        let mut size = vec![0usize; groups];
        for (h, &g) in group.iter().enumerate() {
            if g >= groups {
                return Err(Error(format!(
                    "reference haplotype {h} is in group {g}, but there are {groups} groups"
                )));
            }
            size[g] += 1;
        }
        if let Some(g) = size.iter().position(|&n| n == 0) {
            return Err(Error(format!("group {g} has no reference haplotype")));
        }
        Ok(RefPanel {
            haplotypes: group.len(),
            words: words(group.len()),
            bits: Vec::new(),
            group,
            size,
        })
    }

    /// The panel of `haplotypes` (every site it holds); `group[p]` is the
    /// group of person `p`, so both haplotypes of a person share it.
    pub fn from_haplotypes(
        haplotypes: &Haplotypes,
        person_group: &[usize],
        groups: usize,
    ) -> Result<RefPanel> {
        if person_group.len() != haplotypes.people().len() {
            return Err(Error(format!(
                "{} groups given for {} reference people",
                person_group.len(),
                haplotypes.people().len()
            )));
        }
        let group = person_group.iter().flat_map(|&g| [g, g]).collect();
        let mut panel = RefPanel::new(group, groups)?;
        for m in 0..haplotypes.sites().len() {
            panel.push_site(haplotypes.row(m));
        }
        Ok(panel)
    }

    /// Adds a site: bit `h` of `row` is the allele of haplotype `h`.
    pub fn push_site(&mut self, row: &[u64]) {
        assert_eq!(row.len(), self.words, "a row of {} words", self.words);
        self.bits.extend_from_slice(row);
        // Bits past the last haplotype are never read, but zero them so
        // equal panels compare equal.
        let extra = self.words * 64 - self.haplotypes;
        if extra > 0 {
            let last = self.bits.len() - 1;
            self.bits[last] &= u64::MAX >> extra;
        }
    }

    pub fn sites(&self) -> usize {
        self.bits.len().checked_div(self.words).unwrap_or(0)
    }

    pub fn haplotypes(&self) -> usize {
        self.haplotypes
    }

    pub fn groups(&self) -> usize {
        self.size.len()
    }

    /// The group of haplotype `h`.
    pub fn group_of(&self, h: usize) -> usize {
        self.group[h]
    }

    /// The group of every haplotype.
    pub fn group_all(&self) -> &[usize] {
        &self.group
    }

    /// The haplotypes of each group.
    pub fn sizes(&self) -> &[usize] {
        &self.size
    }

    /// The alleles of every haplotype at site `m`.
    pub fn row(&self, m: usize) -> &[u64] {
        &self.bits[m * self.words..(m + 1) * self.words]
    }

    /// The allele (0 or 1) of haplotype `h` at site `m`.
    pub fn allele(&self, m: usize, h: usize) -> u8 {
        bit(self.row(m), h)
    }
}

/// One haplotype to paint: one allele bit per site.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Bits {
    len: usize,
    words: Vec<u64>,
}

impl Bits {
    /// From alleles, each 0 or 1 (anything else counts as 1).
    pub fn from_alleles(alleles: &[u8]) -> Bits {
        let mut b = Bits::default();
        for &a in alleles {
            b.push(a != 0);
        }
        b
    }

    pub fn push(&mut self, one: bool) {
        if self.len.is_multiple_of(64) {
            self.words.push(0);
        }
        if one {
            self.words[self.len / 64] |= 1 << (self.len % 64);
        }
        self.len += 1;
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The allele at site `m`.
    pub fn get(&self, m: usize) -> u8 {
        assert!(m < self.len, "site {m} of {}", self.len);
        bit(&self.words, m)
    }
}
