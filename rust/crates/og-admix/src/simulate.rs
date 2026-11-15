//! An admixed diploid genome with known truth, the input every engine's
//! truth tests start from (`og ancestry-fine simulate`).
//!
//! Each of the two haplotypes is built on its own along the genetic map:
//!
//! - **Ancestry tracts.** Ancestry events come as a Poisson process of
//!   rate T per Morgan (T is about the generations since admixture). At
//!   each event the new ancestry is drawn from the proportions mu, so it
//!   may be the ancestry the haplotype already has, as in the model the
//!   painting fits. The raw tract lengths are exponential with mean 1/T
//!   Morgan.
//! - **Copying inside a tract.** The haplotype copies one haplotype of a
//!   held-out person of that ancestry, drawn at the start of the tract and
//!   drawn again at copying switches, which come at `switch_rate` per
//!   Morgan (1 by default: one crossover per Morgan per meiosis). A draw
//!   never takes the haplotype the other copy is copying at that moment
//!   when its group has another one, so a group of one held-out person
//!   gives back that person's own genotype instead of a run of
//!   homozygosity.
//!
//! The genome is written as an unphased VCF with sorted alleles (`0/1`,
//! never `1/0`), so nothing in it says which copy an allele came from and
//! Beagle phases it as it would a real person. Every biallelic SNV of the
//! panel is written, the first one at a position where there are two. The
//! truth is a BED of each haplotype's ancestry: per chromosome and
//! haplotype the rows run from the first written site to the last with no
//! gap and no overlap, and two touching rows never share an ancestry. A
//! breakpoint between two sites sits where the map puts the event, kept
//! between the two sites, so every site's truth is the row that holds it.
//!
//! Each chromosome and haplotype draws from its own [`Rng::stream`] of the
//! seed, so a chromosome's genome does not change when other chromosomes
//! are added or left out, and the same seed gives the same bytes.

use crate::map::GeneticMap;
use crate::panel::{PanelReader, bit, words};
use crate::rng::Rng;
use crate::{Error, Result, repr};
use std::io::Write;
use std::path::Path;

/// How a genome is drawn.
#[derive(Debug, Clone, PartialEq)]
pub struct Config {
    /// T, ancestry events per Morgan.
    pub generations: f64,
    /// Copying switches per Morgan inside a tract.
    pub switch_rate: f64,
    /// The ancestries and their proportions mu, in the order the truth
    /// numbers them.
    pub mu: Vec<(String, f64)>,
    pub seed: u64,
}

/// The most events per Morgan T or the switch rate may ask for. A higher
/// rate gives tracts shorter than the gap between two panel sites, and at
/// about 1e16 an event's place in cM stops moving, so the draw never ends.
pub const MAX_PER_MORGAN: f64 = 10_000.0;

impl Config {
    /// Checks the numbers: T above 0, a switch rate of 0 or more, both at
    /// most [`MAX_PER_MORGAN`], and proportions above 0 that sum to 1
    /// (within 1e-6), one per ancestry.
    pub fn check(&self) -> Result<()> {
        if !(self.generations.is_finite() && self.generations > 0.0) {
            return Err(Error(format!(
                "T is {}: ancestry events per Morgan must be above 0",
                self.generations
            )));
        }
        if self.generations > MAX_PER_MORGAN {
            return Err(Error(format!(
                "T is {}: ancestry events per Morgan must be at most {MAX_PER_MORGAN}",
                self.generations
            )));
        }
        if !(self.switch_rate.is_finite() && self.switch_rate >= 0.0) {
            return Err(Error(format!(
                "the switch rate is {}: copying switches per Morgan must be 0 or more",
                self.switch_rate
            )));
        }
        if self.switch_rate > MAX_PER_MORGAN {
            return Err(Error(format!(
                "the switch rate is {}: copying switches per Morgan must be at most {MAX_PER_MORGAN}",
                self.switch_rate
            )));
        }
        if self.mu.is_empty() {
            return Err(Error("no ancestry proportions".into()));
        }
        for (i, (g, w)) in self.mu.iter().enumerate() {
            if self.mu[..i].iter().any(|(h, _)| h == g) {
                return Err(Error(format!("ancestry {} is given twice", repr(g))));
            }
            if !(w.is_finite() && *w > 0.0) {
                return Err(Error(format!(
                    "the proportion of {} is {w}: each must be above 0",
                    repr(g)
                )));
            }
        }
        let total: f64 = self.mu.iter().map(|(_, w)| w).sum();
        if (total - 1.0).abs() > 1e-6 {
            return Err(Error(format!(
                "the ancestry proportions sum to {total}, not 1"
            )));
        }
        Ok(())
    }

    /// Running totals of mu, for [`Rng::pick`].
    fn cumulative(&self) -> Vec<f64> {
        self.mu
            .iter()
            .scan(0.0, |s, (_, w)| {
                *s += w;
                Some(*s)
            })
            .collect()
    }
}

/// The ancestry events of one haplotype along the map: a Poisson process
/// of rate T per Morgan, each event drawing the ancestry from mu.
#[derive(Debug, Clone)]
pub struct AncestryProcess {
    per_cm: f64,
    cumulative: Vec<f64>,
    next_cm: f64,
    ancestry: usize,
}

impl AncestryProcess {
    /// Starts at `start_cm`: draws the first ancestry and the first event.
    /// `cumulative` holds the running totals of mu.
    pub fn new(rng: &mut Rng, generations: f64, cumulative: &[f64], start_cm: f64) -> Self {
        let per_cm = generations / 100.0;
        let ancestry = rng.pick(cumulative);
        AncestryProcess {
            per_cm,
            cumulative: cumulative.to_vec(),
            next_cm: start_cm + rng.exponential(per_cm),
            ancestry,
        }
    }

    pub fn ancestry(&self) -> usize {
        self.ancestry
    }

    /// Where the next event is, in cM.
    pub fn next_cm(&self) -> f64 {
        self.next_cm
    }

    /// Moves to the next event: its place in cM and the ancestry drawn
    /// there.
    pub fn step(&mut self, rng: &mut Rng) -> (f64, usize) {
        let at = self.next_cm;
        self.ancestry = rng.pick(&self.cumulative);
        self.next_cm = at + rng.exponential(self.per_cm);
        (at, self.ancestry)
    }
}

/// One row of the truth: `[start, end)` in 0-based BED coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TruthRow {
    pub start: i64,
    pub end: i64,
    /// 0 or 1.
    pub haplotype: u8,
    /// An index into [`Config::mu`].
    pub ancestry: usize,
}

/// What one chromosome gave.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: String,
    /// Sites written.
    pub sites: u64,
    /// Records of the panel left out: not a biallelic SNV, or a second SNV
    /// at a position.
    pub skipped: u64,
    /// The truth, by start, haplotype 0 before 1.
    pub rows: Vec<TruthRow>,
}

/// One haplotype being built.
struct Copy {
    rng: Rng,
    tracts: AncestryProcess,
    switch_per_cm: f64,
    next_switch: f64,
    /// An index into the haplotypes of the current ancestry's group.
    source: usize,
    start: i64,
    first: usize,
    /// (breakpoint, ancestry from there on).
    changes: Vec<(i64, usize)>,
}

/// A source among `n`, never `exclude` when there is another.
fn draw_source(rng: &mut Rng, n: usize, exclude: Option<usize>) -> usize {
    match exclude {
        Some(x) if n > 1 => {
            let i = rng.below(n as u64 - 1) as usize;
            if i >= x { i + 1 } else { i }
        }
        _ => rng.below(n as u64) as usize,
    }
}

impl Copy {
    fn new(
        mut rng: Rng,
        cfg: &Config,
        cumulative: &[f64],
        groups: &[Vec<usize>],
        start_cm: f64,
        start: i64,
        other: Option<(usize, usize)>,
    ) -> Copy {
        let tracts = AncestryProcess::new(&mut rng, cfg.generations, cumulative, start_cm);
        let a = tracts.ancestry();
        let exclude = other.and_then(|(oa, os)| (oa == a).then_some(os));
        let source = draw_source(&mut rng, groups[a].len(), exclude);
        let switch_per_cm = cfg.switch_rate / 100.0;
        let next_switch = start_cm + rng.exponential(switch_per_cm);
        Copy {
            rng,
            tracts,
            switch_per_cm,
            next_switch,
            source,
            start,
            first: a,
            changes: Vec::new(),
        }
    }

    fn state(&self) -> (usize, usize) {
        (self.tracts.ancestry(), self.source)
    }

    /// Takes every event up to `cm`, the place of the site at `pos`; the
    /// site before it is at `prev`.
    fn advance(
        &mut self,
        cm: f64,
        prev: i64,
        pos: i64,
        map: &GeneticMap,
        groups: &[Vec<usize>],
        other: (usize, usize),
    ) {
        loop {
            let tract_next = self.tracts.next_cm();
            if tract_next.min(self.next_switch) > cm {
                return;
            }
            let at = if tract_next <= self.next_switch {
                let (at, a) = self.tracts.step(&mut self.rng);
                // The site before keeps the old ancestry, this one the new.
                let b = (map.pos_at(at).round() as i64).clamp(prev, pos - 1);
                self.changes.push((b, a));
                at
            } else {
                self.next_switch
            };
            let a = self.tracts.ancestry();
            let exclude = (other.0 == a).then_some(other.1);
            self.source = draw_source(&mut self.rng, groups[a].len(), exclude);
            self.next_switch = at + self.rng.exponential(self.switch_per_cm);
        }
    }

    /// The truth rows up to `end`, empty rows dropped and touching rows of
    /// one ancestry joined.
    fn rows(&self, haplotype: u8, end: i64) -> Vec<TruthRow> {
        let mut out: Vec<TruthRow> = Vec::new();
        let (mut start, mut ancestry) = (self.start, self.first);
        let mut push = |start: i64, end: i64, ancestry: usize| {
            if end <= start {
                return;
            }
            match out.last_mut() {
                Some(last) if last.ancestry == ancestry => last.end = end,
                _ => out.push(TruthRow {
                    start,
                    end,
                    haplotype,
                    ancestry,
                }),
            }
        };
        for &(b, a) in &self.changes {
            if b > start {
                push(start, b, ancestry);
                start = b;
            }
            ancestry = a;
        }
        push(start, end, ancestry);
        out
    }
}

/// Writes the VCF header: one `##contig` per chromosome, the settings, and
/// one sample.
pub fn write_header(
    out: &mut dyn Write,
    sample: &str,
    chroms: &[&str],
    cfg: &Config,
) -> Result<()> {
    let mu: Vec<String> = cfg.mu.iter().map(|(g, w)| format!("{g}:{w}")).collect();
    let mut text = String::from("##fileformat=VCFv4.2\n##source=og ancestry-fine simulate\n");
    text.push_str(&format!(
        "##og_simulate=seed={};generations={};switch_rate={};mu={}\n",
        cfg.seed,
        cfg.generations,
        cfg.switch_rate,
        mu.join(",")
    ));
    for c in chroms {
        text.push_str(&format!("##contig=<ID={c}>\n"));
    }
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    text.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ));
    out.write_all(text.as_bytes())
        .map_err(|e| Error(format!("cannot write the VCF: {e}")))
}

/// Simulates one chromosome from the panel file `panel` and its `map`:
/// writes its records to `vcf` and gives back its truth. `people` holds the
/// held-out people as (sample, ancestry); only those of an ancestry in mu
/// are read, and each ancestry in mu needs at least one.
pub fn simulate_chrom(
    panel: &Path,
    map: &GeneticMap,
    people: &[(String, String)],
    cfg: &Config,
    vcf: &mut dyn Write,
) -> Result<Chrom> {
    cfg.check()?;
    if map.is_empty() {
        return Err(Error(format!(
            "no genetic map markers for {}",
            panel.display()
        )));
    }
    let cumulative = cfg.cumulative();
    let mut ids: Vec<&str> = Vec::new();
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); cfg.mu.len()];
    for (id, group) in people {
        if let Some(g) = cfg.mu.iter().position(|(m, _)| m == group) {
            groups[g].extend([2 * ids.len(), 2 * ids.len() + 1]);
            ids.push(id);
        }
    }
    if let Some(((g, _), _)) = cfg.mu.iter().zip(&groups).find(|(_, h)| h.is_empty()) {
        return Err(Error(format!(
            "no held-out person of ancestry {} to copy",
            repr(g)
        )));
    }
    let mut reader = PanelReader::open(panel, &ids)?;
    let mut row = vec![0u64; words(reader.haplotypes())];
    let mut copies: Vec<Copy> = Vec::with_capacity(2);
    let (mut prev, mut sites, mut doubles) = (0i64, 0u64, 0u64);
    let mut line = Vec::with_capacity(64);
    let mut name = String::new();
    while let Some(site) = reader.next_site(&mut row)? {
        let pos = i64::from(site.pos);
        if sites > 0 && pos == prev {
            doubles += 1;
            continue;
        }
        let cm = map.cm_at(pos);
        if copies.is_empty() {
            name = reader.chrom().unwrap_or("").to_owned();
            for h in 0..2u8 {
                let rng = Rng::stream(cfg.seed, &format!("{name}/{h}"));
                let other = copies.first().map(Copy::state);
                copies.push(Copy::new(
                    rng,
                    cfg,
                    &cumulative,
                    &groups,
                    cm,
                    pos - 1,
                    other,
                ));
            }
        } else {
            for h in 0..2 {
                let other = copies[1 - h].state();
                copies[h].advance(cm, prev, pos, map, &groups, other);
            }
        }
        let mut alleles = [0u8; 2];
        for (h, copy) in copies.iter().enumerate() {
            let (a, s) = copy.state();
            alleles[h] = bit(&row, groups[a][s]);
        }
        let gt: &[u8] = match alleles[0] + alleles[1] {
            0 => b"0/0",
            1 => b"0/1",
            _ => b"1/1",
        };
        line.clear();
        line.extend_from_slice(name.as_bytes());
        line.extend_from_slice(format!("\t{pos}\t.\t").as_bytes());
        line.extend_from_slice(&[site.ref_, b'\t', site.alt]);
        line.extend_from_slice(b"\t.\tPASS\t.\tGT\t");
        line.extend_from_slice(gt);
        line.push(b'\n');
        vcf.write_all(&line)
            .map_err(|e| Error(format!("cannot write the VCF: {e}")))?;
        prev = pos;
        sites += 1;
    }
    if sites == 0 {
        return Err(Error(format!(
            "{}: no biallelic SNV to copy",
            panel.display()
        )));
    }
    let mut rows: Vec<TruthRow> = copies
        .iter()
        .enumerate()
        .flat_map(|(h, c)| c.rows(h as u8, prev))
        .collect();
    rows.sort_by_key(|r| (r.start, r.haplotype));
    Ok(Chrom {
        name,
        sites,
        skipped: reader.skipped() + doubles,
        rows,
    })
}

/// Writes the truth BED of every chromosome: a header line, then
/// `chrom start end haplotype ancestry` rows.
pub fn write_truth(out: &mut dyn Write, chroms: &[Chrom], cfg: &Config) -> Result<()> {
    let mut text = String::from("#chrom\tstart\tend\thaplotype\tancestry\n");
    for c in chroms {
        for r in &c.rows {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                c.name, r.start, r.end, r.haplotype, cfg.mu[r.ancestry].0
            ));
        }
    }
    out.write_all(text.as_bytes())
        .map_err(|e| Error(format!("cannot write the truth BED: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    /// Whether `lengths` (cM) follow the exponential law of mean 100/T cM:
    /// the mean within 5%, and the share longer than the mean within 0.02
    /// of 1/e.
    fn exponential_law(lengths: &[f64], generations: f64) -> std::result::Result<(), String> {
        let mean_want = 100.0 / generations;
        let n = lengths.len() as f64;
        let mean = lengths.iter().sum::<f64>() / n;
        let above = lengths.iter().filter(|&&l| l > mean_want).count() as f64 / n;
        if (mean / mean_want - 1.0).abs() > 0.05 {
            return Err(format!("mean tract {mean} cM, the law says {mean_want}"));
        }
        if (above - (-1.0f64).exp()).abs() > 0.02 {
            return Err(format!(
                "{above} of tracts above the mean, the law says 1/e"
            ));
        }
        Ok(())
    }

    fn tract_lengths(generations: f64, seed: u64, n: usize) -> Vec<f64> {
        let mut rng = Rng::new(seed);
        let mut p = AncestryProcess::new(&mut rng, generations, &[0.3, 1.0], 0.0);
        let mut last = 0.0;
        (0..n)
            .map(|_| {
                let (at, _) = p.step(&mut rng);
                let l = at - last;
                last = at;
                l
            })
            .collect()
    }

    #[test]
    fn tract_lengths_follow_the_exponential_law() {
        let t = 10.0;
        for seed in 1..=5 {
            let drawn_t = if crate::PLANTED { 2.0 * t } else { t };
            let lengths = tract_lengths(drawn_t, seed, 10_000);
            exponential_law(&lengths, t).unwrap();
        }
    }

    #[test]
    fn the_law_check_fails_when_t_is_doubled() {
        for seed in 1..=5 {
            let lengths = tract_lengths(20.0, seed, 10_000);
            let err = exponential_law(&lengths, 10.0).unwrap_err();
            assert!(err.starts_with("mean tract"), "{err}");
        }
    }

    /// A panel where group A carries REF and group B carries ALT at every
    /// site, so a site's ALT count is the number of copies of B there; plus
    /// a person of no group. `n` sites, 1 kb apart, from position 1001.
    fn two_group_panel(d: &Dir, n: u32) -> std::path::PathBuf {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\tB1\tX\tA2\tB2\n",
        );
        for i in 1..=n {
            text.push_str(&format!(
                "chr7\t{}\t.\tC\tT\t.\t.\t.\tGT\t0|0\t1|1\t0|1\t0|0\t1|1\n",
                1000 + i * 1000
            ));
            if i % 97 == 0 {
                // A second SNV at the position, left out.
                text.push_str(&format!(
                    "chr7\t{}\t.\tC\tG\t.\t.\t.\tGT\t0|0\t1|1\t0|1\t0|0\t1|1\n",
                    1000 + i * 1000
                ));
            }
        }
        d.put("chr7.vcf.gz", &text)
    }

    /// 100 cM per Mb over the panel: 3 Morgans for 3,000 sites.
    fn steep_map() -> GeneticMap {
        GeneticMap::from_markers(vec![(1, 0.0), (3_001_001, 300.0)]).unwrap()
    }

    fn people() -> Vec<(String, String)> {
        [
            ("A1", "A"),
            ("B1", "B"),
            ("X", "Z"),
            ("A2", "A"),
            ("B2", "B"),
        ]
        .iter()
        .map(|(p, g)| ((*p).to_owned(), (*g).to_owned()))
        .collect()
    }

    fn cfg(seed: u64) -> Config {
        Config {
            generations: 5.0,
            switch_rate: 1.0,
            mu: vec![("A".into(), 0.6), ("B".into(), 0.4)],
            seed,
        }
    }

    /// Every haplotype's rows run from the first site to the last with no
    /// gap or overlap, and touching rows differ in ancestry.
    fn check_cover(rows: &[TruthRow], start: i64, end: i64) -> std::result::Result<(), String> {
        for h in 0..2u8 {
            let mine: Vec<&TruthRow> = rows.iter().filter(|r| r.haplotype == h).collect();
            let mut at = start;
            let mut last = None;
            for r in &mine {
                if r.start != at {
                    return Err(format!(
                        "haplotype {h}: a row starts at {}, not {at}",
                        r.start
                    ));
                }
                if r.end <= r.start {
                    return Err(format!("haplotype {h}: an empty row at {}", r.start));
                }
                if last == Some(r.ancestry) {
                    return Err(format!(
                        "haplotype {h}: two touching rows of one ancestry at {at}"
                    ));
                }
                at = r.end;
                last = Some(r.ancestry);
            }
            if at != end {
                return Err(format!("haplotype {h}: the rows end at {at}, not {end}"));
            }
        }
        Ok(())
    }

    fn ancestry_at(rows: &[TruthRow], h: u8, pos: i64) -> usize {
        rows.iter()
            .find(|r| r.haplotype == h && r.start < pos && pos <= r.end)
            .map(|r| r.ancestry)
            .expect("a row holds every site")
    }

    #[test]
    fn truth_covers_the_chromosome_and_explains_every_genotype() {
        let d = Dir::new("sim-truth");
        let panel = two_group_panel(&d, 3000);
        let map = steep_map();
        let (mut bp_a, mut bp_all) = (0i64, 0i64);
        let mut breaks = 0usize;
        for seed in 1..=20 {
            let mut vcf = Vec::new();
            let chrom = simulate_chrom(&panel, &map, &people(), &cfg(seed), &mut vcf).unwrap();
            assert_eq!(
                (chrom.name.as_str(), chrom.sites, chrom.skipped),
                ("chr7", 3000, 30)
            );
            let mut rows = chrom.rows.clone();
            if crate::PLANTED && rows.len() > 2 {
                rows.remove(1);
            }
            check_cover(&rows, 2000 - 1, 3_001_000).unwrap();
            for line in String::from_utf8(vcf).unwrap().lines() {
                let c: Vec<&str> = line.split('\t').collect();
                let pos: i64 = c[1].parse().unwrap();
                let alt = c[9].bytes().filter(|&b| b == b'1').count();
                let want = (0..2).filter(|&h| ancestry_at(&rows, h, pos) == 1).count();
                assert_eq!(alt, want, "seed {seed}, {line}");
                assert!(["0/0", "0/1", "1/1"].contains(&c[9]), "{line}");
                assert_eq!(c[3..5], ["C", "T"]);
            }
            for r in &rows {
                bp_all += r.end - r.start;
                if r.ancestry == 0 {
                    bp_a += r.end - r.start;
                }
            }
            breaks += rows.len() - 2;
        }
        // 40 haplotypes of 3 Morgans: share of A near mu, and the joined
        // breakpoints near T (1 - sum mu^2) = 2.4 per Morgan.
        let share = bp_a as f64 / bp_all as f64;
        assert!((share - 0.6).abs() < 0.06, "share of A {share}");
        let per_morgan = breaks as f64 / 120.0;
        assert!(
            (per_morgan - 2.4).abs() < 0.5,
            "{per_morgan} breakpoints per Morgan"
        );
    }

    #[test]
    fn the_cover_check_finds_a_gap_an_overlap_and_a_split() {
        let r = |start, end, haplotype, ancestry| TruthRow {
            start,
            end,
            haplotype,
            ancestry,
        };
        let good = [r(0, 5, 0, 0), r(0, 10, 1, 1), r(5, 10, 0, 1)];
        check_cover(&good, 0, 10).unwrap();
        let gap = [r(0, 4, 0, 0), r(0, 10, 1, 1), r(5, 10, 0, 1)];
        assert!(
            check_cover(&gap, 0, 10)
                .unwrap_err()
                .contains("starts at 5, not 4")
        );
        let overlap = [r(0, 6, 0, 0), r(0, 10, 1, 1), r(5, 10, 0, 1)];
        assert!(check_cover(&overlap, 0, 10).is_err());
        let split = [r(0, 5, 0, 1), r(0, 10, 1, 1), r(5, 10, 0, 1)];
        assert!(
            check_cover(&split, 0, 10)
                .unwrap_err()
                .contains("one ancestry")
        );
        assert!(check_cover(&good, 0, 11).unwrap_err().contains("not 11"));
    }

    #[test]
    fn one_held_out_person_comes_back_as_themselves() {
        let d = Dir::new("sim-self");
        let mut text = String::from(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tW\n",
        );
        let gts = ["0|1", "1|0", "1|1", "0|0", "0|1"];
        for i in 0..2000u32 {
            text.push_str(&format!(
                "chr3\t{}\t.\tA\tG\t.\t.\t.\tGT\t{}\n",
                100 + i * 1000,
                gts[(i * 7 % 5) as usize]
            ));
        }
        let panel = d.put("chr3.vcf.gz", &text);
        let map = GeneticMap::from_markers(vec![(1, 0.0), (2_000_100, 200.0)]).unwrap();
        let cfg = Config {
            generations: 6.0,
            switch_rate: 5.0,
            mu: vec![("WEA".into(), 1.0)],
            seed: 9,
        };
        let mut vcf = Vec::new();
        let chrom =
            simulate_chrom(&panel, &map, &[("W".into(), "WEA".into())], &cfg, &mut vcf).unwrap();
        assert_eq!(chrom.rows.len(), 2, "one row per haplotype");
        for (i, line) in String::from_utf8(vcf).unwrap().lines().enumerate() {
            let want = match gts[i * 7 % 5] {
                "0|0" => "0/0",
                "1|1" => "1/1",
                _ => "0/1",
            };
            assert!(line.ends_with(want), "{line}");
        }
    }

    fn whole_run(seed: u64) -> (Vec<u8>, Vec<u8>) {
        let d = Dir::new(&format!("sim-run-{seed}"));
        let panel = two_group_panel(&d, 500);
        let map = steep_map();
        let cfg = cfg(seed);
        let mut vcf = Vec::new();
        write_header(&mut vcf, "SIM", &["chr7"], &cfg).unwrap();
        let chrom = simulate_chrom(&panel, &map, &people(), &cfg, &mut vcf).unwrap();
        let mut bed = Vec::new();
        write_truth(&mut bed, &[chrom], &cfg).unwrap();
        (vcf, bed)
    }

    #[test]
    fn the_same_seed_gives_the_same_bytes() {
        let a = whole_run(4);
        let b = whole_run(if crate::PLANTED { 5 } else { 4 });
        assert_eq!(a, b);
        let c = whole_run(5);
        assert_ne!(a.1, c.1, "another seed, another truth");
        let vcf = String::from_utf8(a.0).unwrap();
        assert!(vcf.starts_with("##fileformat=VCFv4.2\n##source=og ancestry-fine simulate\n##og_simulate=seed=4;generations=5;switch_rate=1;mu=A:0.6,B:0.4\n##contig=<ID=chr7>\n"));
        assert!(vcf.contains("\nchr7\t2000\t.\tC\tT\t.\tPASS\t.\tGT\t"));
        let bed = String::from_utf8(a.1).unwrap();
        assert!(
            bed.starts_with("#chrom\tstart\tend\thaplotype\tancestry\nchr7\t1999\t"),
            "{bed}"
        );
    }

    #[test]
    fn settings_are_checked() {
        let d = Dir::new("sim-check");
        let panel = two_group_panel(&d, 10);
        let map = steep_map();
        let run = |cfg: Config, people: &[(String, String)]| {
            simulate_chrom(&panel, &map, people, &cfg, &mut Vec::new())
                .unwrap_err()
                .0
        };
        let mut c = cfg(1);
        c.generations = 0.0;
        assert!(run(c, &people()).contains("must be above 0"));
        let mut c = cfg(1);
        c.generations = 1e300;
        assert!(run(c, &people()).contains("at most 10000"));
        let mut c = cfg(1);
        c.switch_rate = -1.0;
        assert!(run(c, &people()).contains("must be 0 or more"));
        let mut c = cfg(1);
        c.switch_rate = 20_000.0;
        assert!(run(c, &people()).contains("at most 10000"));
        let mut c = cfg(1);
        c.mu[1].1 = 0.5;
        assert!(run(c, &people()).contains("sum to 1.1, not 1"));
        let mut c = cfg(1);
        c.mu[1].0 = "A".into();
        assert!(run(c, &people()).contains("ancestry 'A' is given twice"));
        let mut c = cfg(1);
        c.mu = vec![("A".into(), 0.5), ("Q".into(), 0.5)];
        assert!(run(c, &people()).contains("no held-out person of ancestry 'Q'"));
        let e = simulate_chrom(
            &panel,
            &GeneticMap::default(),
            &people(),
            &cfg(1),
            &mut Vec::new(),
        )
        .unwrap_err()
        .0;
        assert!(e.starts_with("no genetic map markers"), "{e}");
    }
}
