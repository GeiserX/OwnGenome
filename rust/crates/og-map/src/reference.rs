//! GRCh38's chromosomes from UCSC's band file (`assets/reference/cytoBand.hg38.txt`):
//! the bands, each chromosome's length (its last band's end), its centromere
//! (the `acen` bands) and its masked spans (the `acen`, `gvar` and `stalk`
//! bands, neighbours merged). One file, so the lengths and the masked spans
//! cannot drift from the bands the map draws.

use std::fmt;
use std::io::BufRead;
use std::path::Path;

/// The chromosomes the map draws, in display order.
pub const CHROMS: [&str; 24] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY",
];

/// The band file under `--assets`.
pub const BAND_FILE: &str = "reference/cytoBand.hg38.txt";

/// One band: 0-based start, exclusive end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Band {
    pub s: u64,
    pub e: u64,
    pub name: String,
    pub stain: String,
}

/// One chromosome of the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: &'static str,
    pub len: u64,
    /// The span of its `acen` bands.
    pub cen: (u64, u64),
    /// `acen`, `gvar` and `stalk` bands, touching neighbours merged.
    pub masked: Vec<(u64, u64)>,
    /// Sorted by start.
    pub bands: Vec<Band>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// In [`CHROMS`] order, all 24.
    pub chroms: Vec<Chrom>,
}

/// Why the band file could not be read: `file:line: reason`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefError(pub String);

impl fmt::Display for RefError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RefError {}

/// The index of `name` in [`CHROMS`].
pub fn chrom_index(name: &str) -> Option<usize> {
    CHROMS.iter().position(|&c| c == name)
}

fn masked_stain(stain: &str) -> bool {
    matches!(stain, "acen" | "gvar" | "stalk")
}

impl Reference {
    /// Reads `<assets>/reference/cytoBand.hg38.txt`.
    pub fn load(assets: &Path) -> Result<Reference, RefError> {
        let path = assets.join(BAND_FILE);
        let f =
            std::fs::File::open(&path).map_err(|e| RefError(format!("{}: {e}", path.display())))?;
        Reference::parse(std::io::BufReader::new(f), "cytoBand.hg38.txt")
    }

    /// Parses band rows (chrom, start, end, name, stain); rows on other
    /// contigs are skipped. `name` labels the errors.
    pub fn parse(input: impl BufRead, name: &str) -> Result<Reference, RefError> {
        let mut bands: Vec<Vec<Band>> = vec![Vec::new(); CHROMS.len()];
        for (i, line) in input.lines().enumerate() {
            let at = |why: String| RefError(format!("{name}:{}: {why}", i + 1));
            let line = line.map_err(|e| at(e.to_string()))?;
            let line = line.trim_end_matches('\r');
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 5 {
                return Err(at(format!("expected 5 columns, found {}", f.len())));
            }
            let Some(k) = chrom_index(f[0]) else { continue };
            let num = |s: &str| {
                s.parse::<u64>()
                    .map_err(|_| at(format!("not a position: {s:?}")))
            };
            let (s, e) = (num(f[1])?, num(f[2])?);
            if e <= s {
                return Err(at(format!("end {e} is not after start {s}")));
            }
            bands[k].push(Band {
                s,
                e,
                name: f[3].to_owned(),
                stain: f[4].to_owned(),
            });
        }
        let mut chroms = Vec::with_capacity(CHROMS.len());
        for (k, mut b) in bands.into_iter().enumerate() {
            let cname = CHROMS[k];
            if b.is_empty() {
                return Err(RefError(format!("{name}: no band on {cname}")));
            }
            b.sort_by_key(|x| (x.s, x.e));
            let len = b.iter().map(|x| x.e).max().unwrap_or(0);
            let acen: Vec<&Band> = b.iter().filter(|x| x.stain == "acen").collect();
            let cen = match (acen.first(), acen.last()) {
                (Some(a), Some(z)) => (a.s, z.e),
                _ => return Err(RefError(format!("{name}: no centromere band on {cname}"))),
            };
            let mut masked: Vec<(u64, u64)> = Vec::new();
            for x in b.iter().filter(|x| masked_stain(&x.stain)) {
                match masked.last_mut() {
                    Some(last) if x.s <= last.1 => last.1 = last.1.max(x.e),
                    _ => masked.push((x.s, x.e)),
                }
            }
            chroms.push(Chrom {
                name: cname,
                len,
                cen,
                masked,
                bands: b,
            });
        }
        Ok(Reference { chroms })
    }

    /// The chromosome named `name`.
    pub fn chrom(&self, name: &str) -> Option<&Chrom> {
        chrom_index(name).map(|k| &self.chroms[k])
    }
}

impl Chrom {
    /// How many bases of `[s, e)` lie in the masked spans.
    pub fn masked_overlap(&self, s: u64, e: u64) -> u64 {
        self.masked
            .iter()
            .map(|&(a, b)| e.min(b).saturating_sub(s.max(a)))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn asset() -> Reference {
        let assets = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
        Reference::load(&assets).expect("the band asset reads")
    }

    /// GRCh38 primary assembly lengths, copied from bin/report_charts.py:56-61
    /// (`GRCH38`, UCSC hg38 chrom.sizes).
    const GRCH38: [(&str, u64); 24] = [
        ("chr1", 248956422),
        ("chr2", 242193529),
        ("chr3", 198295559),
        ("chr4", 190214555),
        ("chr5", 181538259),
        ("chr6", 170805979),
        ("chr7", 159345973),
        ("chr8", 145138636),
        ("chr9", 138394717),
        ("chr10", 133797422),
        ("chr11", 135086622),
        ("chr12", 133275309),
        ("chr13", 114364328),
        ("chr14", 107043718),
        ("chr15", 101991189),
        ("chr16", 90338345),
        ("chr17", 83257441),
        ("chr18", 80373285),
        ("chr19", 58617616),
        ("chr20", 64444167),
        ("chr21", 46709983),
        ("chr22", 50818468),
        ("chrX", 156040895),
        ("chrY", 57227415),
    ];

    /// MASKED_BANDS, copied from bin/collect_summary.py:2477-2483.
    const MASKED_BANDS: &str = concat!(
        "chr1:121700000-143200000 chr2:91800000-96000000 chr3:87800000-98600000 chr4:48200000-51800000 ",
        "chr5:46100000-51400000 chr6:58500000-62600000 chr7:58100000-62100000 chr8:43200000-47200000 ",
        "chr9:42200000-61500000 chr10:38000000-41600000 chr11:51000000-55800000 chr12:33200000-37800000 ",
        "chr13:0-18900000 chr14:0-18200000 chr15:0-20500000 chr16:35300000-47000000 ",
        "chr17:22700000-27400000 chr18:15400000-21500000 chr19:19900000-31900000 chr20:25700000-30400000 ",
        "chr21:0-13000000 chr22:0-17400000 chrX:58100000-63800000 chrY:10300000-10600000 ",
        "chrY:26600000-57227415"
    );

    #[test]
    fn lengths_equal_report_charts() {
        let r = asset();
        let got: Vec<(&str, u64)> = r.chroms.iter().map(|c| (c.name, c.len)).collect();
        assert_eq!(got, GRCH38);
        assert_eq!(r.chroms.iter().map(|c| c.bands.len()).sum::<usize>(), 862);
    }

    #[test]
    fn masked_spans_equal_collect_summary() {
        let r = asset();
        let got: Vec<String> = r
            .chroms
            .iter()
            .flat_map(|c| {
                c.masked
                    .iter()
                    .map(move |(s, e)| format!("{}:{s}-{e}", c.name))
            })
            .collect();
        assert_eq!(got.join(" "), MASKED_BANDS);
    }

    #[test]
    fn centromeres_and_overlap() {
        let r = asset();
        let c1 = r.chrom("chr1").unwrap();
        assert_eq!(c1.cen, (121700000, 125100000));
        assert_eq!(c1.masked_overlap(121000000, 122000000), 300000);
        assert_eq!(c1.masked_overlap(0, 1000), 0);
        assert!(r.chrom("chrM").is_none());
    }

    #[test]
    fn bad_rows_name_their_line() {
        let e =
            Reference::parse("chr1\t0\t5\tp1\tgneg\nchr1\t5\n".as_bytes(), "b.txt").unwrap_err();
        assert_eq!(e.0, "b.txt:2: expected 5 columns, found 2");
        let e = Reference::parse("chr1\t9\t5\tp1\tgneg\n".as_bytes(), "b.txt").unwrap_err();
        assert_eq!(e.0, "b.txt:1: end 5 is not after start 9");
    }
}
