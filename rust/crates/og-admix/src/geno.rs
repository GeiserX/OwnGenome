//! `model/ref.geno`: the reference people's genotypes at the model's SNVs,
//! 2 bits each, so a bootstrap can resample people and count their alleles
//! again without the panel.
//!
//! Layout, version 1, all numbers little-endian:
//!
//! | bytes | what |
//! |---|---|
//! | 8 | magic `OGREFGEN` |
//! | 4 | version, `1` |
//! | 8 | J, the SNVs |
//! | 8 | N, the people |
//! | J x ceil(N / 4) | one row per SNV, in `ref.bim` order |
//!
//! In a row, person i is bits `2 (i mod 4)` and up of byte `i / 4`: the
//! copies of the `.bim`'s A1 allele (0, 1 or 2), or 3 for a missing
//! genotype. The unused bits of a row's last byte are 0. About 1,000
//! people at 100,000 SNVs is 25 MB.

use crate::{Error, Result, cannot_read, cannot_write};
use std::io::{self, Read, Write};
use std::path::Path;

pub const MAGIC: [u8; 8] = *b"OGREFGEN";
pub const VERSION: u32 = 1;
/// The code of a missing genotype.
pub const MISSING: u8 = 3;
const HEADER: usize = 8 + 4 + 8 + 8;

/// Genotypes of N people at J SNVs, packed 4 to a byte.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefGeno {
    people: usize,
    sites: usize,
    row: usize,
    data: Vec<u8>,
}

impl RefGeno {
    /// No SNV yet, for `people` people.
    pub fn new(people: usize) -> RefGeno {
        RefGeno {
            people,
            sites: 0,
            row: people.div_ceil(4),
            data: Vec::new(),
        }
    }

    pub fn people(&self) -> usize {
        self.people
    }

    pub fn sites(&self) -> usize {
        self.sites
    }

    /// Adds the next SNV: one genotype per person, the copies of A1, `None`
    /// when missing.
    pub fn push_site(&mut self, genotypes: &[Option<u8>]) -> Result<()> {
        if genotypes.len() != self.people {
            return Err(Error(format!(
                "a ref.geno row holds {} genotypes, not {}",
                self.people,
                genotypes.len()
            )));
        }
        let start = self.data.len();
        self.data.resize(start + self.row, 0);
        for (i, g) in genotypes.iter().enumerate() {
            let code = match *g {
                None => MISSING,
                Some(c @ 0..=2) => c,
                Some(c) => {
                    self.data.truncate(start);
                    return Err(Error(format!(
                        "a genotype counts 0, 1 or 2 copies of A1, not {c}"
                    )));
                }
            };
            self.data[start + i / 4] |= code << (2 * (i % 4));
        }
        self.sites += 1;
        Ok(())
    }

    /// The copies of A1 person `person` carries at SNV `site`, `None` when
    /// missing.
    pub fn get(&self, site: usize, person: usize) -> Option<u8> {
        assert!(
            site < self.sites && person < self.people,
            "outside ref.geno"
        );
        let code = (self.data[site * self.row + person / 4] >> (2 * (person % 4))) & 3;
        (code != MISSING).then_some(code)
    }

    /// At SNV `site`, the copies of A1 among `people` (a person listed twice
    /// counts twice, as a bootstrap draw does) and the people called.
    pub fn a1_count(&self, site: usize, people: &[usize]) -> (u32, u32) {
        people
            .iter()
            .filter_map(|&p| self.get(site, p))
            .fold((0, 0), |(copies, called), g| {
                (copies + u32::from(g), called + 1)
            })
    }

    /// Writes the file's bytes.
    pub fn write(&self, out: &mut impl Write) -> io::Result<()> {
        out.write_all(&MAGIC)?;
        out.write_all(&VERSION.to_le_bytes())?;
        out.write_all(&(self.sites as u64).to_le_bytes())?;
        out.write_all(&(self.people as u64).to_le_bytes())?;
        out.write_all(&self.data)
    }

    /// Reads the file's bytes: the magic, the version, a body of exactly
    /// J rows and nothing after it, and 0 in every unused bit.
    pub fn read(input: &mut impl Read) -> Result<RefGeno> {
        let bad = |what: &str| Error(format!("not a ref.geno file: {what}"));
        let mut head = [0u8; HEADER];
        input
            .read_exact(&mut head)
            .map_err(|_| bad("shorter than its header"))?;
        if head[..8] != MAGIC {
            return Err(bad("no OGREFGEN magic"));
        }
        let word = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().expect("8 bytes"));
        let version = u32::from_le_bytes(head[8..12].try_into().expect("4 bytes"));
        if version != VERSION {
            return Err(Error(format!(
                "ref.geno version {version}: this og reads version {VERSION}; run setup.sh --ancestry-fine-panel again"
            )));
        }
        let (sites, people) = (word(12), word(20));
        let size = usize::try_from(sites)
            .ok()
            .zip(usize::try_from(people).ok())
            .and_then(|(j, n)| j.checked_mul(n.div_ceil(4)).map(|s| (j, n, s)));
        let Some((sites, people, size)) = size else {
            return Err(bad("its size does not fit in memory"));
        };
        let mut data = Vec::new();
        input
            .take((size as u64).saturating_add(1))
            .read_to_end(&mut data)
            .map_err(|e| bad(&e.to_string()))?;
        if data.len() != size {
            return Err(bad(&format!(
                "{} SNVs of {} people take {size} bytes after the header, not {}{}",
                sites,
                people,
                data.len().min(size),
                if data.len() > size { " and more" } else { "" }
            )));
        }
        let geno = RefGeno {
            people,
            sites,
            row: people.div_ceil(4),
            data,
        };
        let used = 2 * (people % 4);
        if used > 0 && geno.row > 0 {
            let mask = !0u8 << used;
            if let Some(j) = (0..sites).find(|j| geno.data[j * geno.row + geno.row - 1] & mask != 0)
            {
                return Err(bad(&format!("unused bits set in the row of SNV {}", j + 1)));
            }
        }
        Ok(geno)
    }

    /// Writes `path` in full under a temporary name, then renames it into
    /// place.
    pub fn save(&self, path: &Path) -> Result<()> {
        let mut f = og_io::AtomicFile::create(path).map_err(|e| cannot_write(path, &e))?;
        self.write(&mut f).map_err(|e| cannot_write(path, &e))?;
        f.commit().map_err(|e| cannot_write(path, &e))
    }

    pub fn load(path: &Path) -> Result<RefGeno> {
        let f = std::fs::File::open(path).map_err(|e| cannot_read(path, &e))?;
        RefGeno::read(&mut io::BufReader::new(f))
            .map_err(|e| Error(format!("{}: {e}", path.display())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;
    use crate::testutil::Dir;

    fn random(rng: &mut Rng, people: usize, sites: usize) -> (RefGeno, Vec<Vec<Option<u8>>>) {
        let mut g = RefGeno::new(people);
        let mut rows = Vec::new();
        for _ in 0..sites {
            let row: Vec<Option<u8>> = (0..people)
                .map(|_| match rng.below(4) {
                    3 => None,
                    c => Some(c as u8),
                })
                .collect();
            g.push_site(&row).unwrap();
            rows.push(row);
        }
        (g, rows)
    }

    #[test]
    fn round_trips_random_matrices() {
        let mut rng = Rng::new(11);
        let d = Dir::new("geno");
        // Every remainder of N mod 4, an empty matrix and a single person.
        for (people, sites) in [(0, 3), (1, 1), (4, 0), (5, 17), (6, 9), (7, 33), (1001, 40)] {
            let (g, rows) = random(&mut rng, people, sites);
            assert_eq!((g.people(), g.sites()), (people, sites));
            let mut bytes = Vec::new();
            g.write(&mut bytes).unwrap();
            assert_eq!(bytes.len(), HEADER + sites * people.div_ceil(4));
            if crate::PLANTED && bytes.len() > HEADER {
                bytes[HEADER] ^= 0b01;
            }
            let back = RefGeno::read(&mut bytes.as_slice()).unwrap();
            assert_eq!(back, g, "{people} people, {sites} SNVs");
            for (j, row) in rows.iter().enumerate() {
                for (i, want) in row.iter().enumerate() {
                    assert_eq!(back.get(j, i), *want);
                }
            }
            let path = d.0.join(format!("ref-{people}-{sites}.geno"));
            g.save(&path).unwrap();
            assert_eq!(RefGeno::load(&path).unwrap(), g);
        }
    }

    #[test]
    fn counts_a1_with_repeats_and_gaps() {
        let mut g = RefGeno::new(3);
        g.push_site(&[Some(2), None, Some(1)]).unwrap();
        assert_eq!(g.a1_count(0, &[0, 1, 2]), (3, 2));
        assert_eq!(g.a1_count(0, &[0, 0, 2]), (5, 3));
        assert!(g.push_site(&[Some(3), None, None]).is_err());
        assert!(g.push_site(&[Some(0)]).is_err());
        assert_eq!(g.sites(), 1, "a refused row adds nothing");
    }

    #[test]
    fn refuses_damaged_files() {
        let mut g = RefGeno::new(5);
        g.push_site(&[Some(0), Some(1), Some(2), None, Some(1)])
            .unwrap();
        g.push_site(&[None; 5]).unwrap();
        let mut good = Vec::new();
        g.write(&mut good).unwrap();
        let read = |b: &[u8]| RefGeno::read(&mut &b[..]).unwrap_err().0;
        assert!(read(&good[..10]).contains("shorter than its header"));
        let mut b = good.clone();
        b[0] = b'X';
        assert!(read(&b).contains("no OGREFGEN magic"));
        let mut b = good.clone();
        b[8] = 2;
        assert!(read(&b).contains("ref.geno version 2"));
        assert!(read(&good[..good.len() - 1]).contains("take 4 bytes after the header, not 3"));
        let mut b = good.clone();
        b.push(0);
        assert!(read(&b).contains("not 4 and more"));
        let mut b = good.clone();
        // Person 5 is in bits 0-1 of the row's second byte; set bit 2.
        b[HEADER + 1] |= 0b100;
        assert!(read(&b).contains("unused bits set in the row of SNV 1"));
        let mut b = good.clone();
        b[12..20].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(read(&b).contains("does not fit in memory"));
        assert_eq!(RefGeno::read(&mut &good[..]).unwrap(), g);
    }
}
