//! Step 44, ancient similarity (`docs/44-ancient.md`): where a genome falls
//! among ancient genomes on a PCA of present-day West Eurasians (smartpca,
//! EIGENSOFT), with the ancient samples of the Allen Ancient DNA Resource
//! (AADR) projected onto the same axes, as `og ancient <verb>`.
//!
//! - `select` (setup): reads the AADR `.anno` and `.ind` files and chooses
//!   the people: the present-day West Eurasians of [`MODERN_GROUPS`] build
//!   the axes, the ancient West Eurasians with enough sites are projected.
//! - `lift` (setup): moves the 1240K sites of the AADR `.snp` file from
//!   GRCh37 to GRCh38 with UCSC's chain file ([`og_intervals::chain`]) and
//!   keeps a site only when one of its two alleles is the GRCh38 base.
//! - `person-geno` (step 44): the person's genotype at each site, appended
//!   as one more column to the present-day panel's EIGENSTRAT `.geno`,
//!   streamed row by row.
//! - `align` (setup, step 44): checks that the axes of two smartpca runs on
//!   the same panel agree, up to the sign of each.
//! - `report` (step 44): the nearest ancient populations, the summary and
//!   the plot.
//!
//! Each command writes what `bin/ancient_report.py` wrote (checked by
//! `tests/parity`). Deliberate changes: every output is written under a
//! temporary name and renamed into place (the Python wrote the population
//! list, `ancients.tsv`, the stats and the text report in place); a gzip
//! input is told by its first bytes, not its name; a `.fai` that does not
//! match its FASTA stops `lift`; and malformed input is an `ERROR:` line
//! and exit 1, never a traceback.

use super::{ArgError, Args, Fail, Pending, Result, Spec, cannot_read, int, repr, split_ws};
use og_intervals::chain::Chain;
use og_intervals::fasta::Indexed;
use og_io::{Lines, Utf8};
use pycompat::fmt::{fixed, int_grouped, int_str, ljust, rjust, round_int};
use pycompat::text::{html_escape, strip, upper};
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Present-day groups of the AADR 1240K file (Group ID, as written in the
/// `.anno` file of v66.p1) that build the axes: Europe, the Caucasus and
/// the Near East. North African groups and groups with recent East
/// Eurasian ancestry are left out, so the axes describe West Eurasian
/// variation only.
pub const MODERN_GROUPS: [&str; 41] = [
    "Abkhasian",
    "Adygei",
    "Albanian",
    "Armenian",
    "Basque",
    "BedouinA",
    "BedouinB",
    "Bulgarian",
    "CEU",
    "Chechen",
    "Czech",
    "Druze",
    "English",
    "Estonian",
    "FIN",
    "Finnish",
    "French",
    "GBR",
    "Georgian",
    "Greek",
    "Greek_Crete",
    "Hungarian",
    "IBS",
    "Icelandic",
    "Iranian",
    "Italian_North",
    "Jew_Iraqi",
    "Jordanian",
    "Lezgin",
    "Norwegian",
    "Orcadian",
    "Ossetian",
    "Palestinian",
    "Polish",
    "Russian",
    "Sardinian",
    "Scottish",
    "Spanish",
    "Spanish_North",
    "TSI",
    "Turkish",
];
const MODERN_DATA_TYPE: &str = "Shotgun.diploid";
/// Ancient samples projected: dated, passing the AADR's own assessment,
/// with this many of the 1240K autosomal sites covered, inside West Eurasia.
const MIN_ANCIENT_SNPS: i64 = 100_000;
const LAT: (f64, f64) = (25.0, 75.0);
const LON: (f64, f64) = (-25.0, 75.0);
const PASSING: [&str; 3] = ["PASS", "PROVISIONAL_PASS", "MERGE_PASS"];
/// The label of a projected ancient person in the `.ind` and `.evec` files.
pub const ANCIENT_LABEL: &str = "Ancient";
/// The label of the person in the run's `.evec`.
pub const PERSON_LABEL: &str = "Person";
/// PCs the distance uses, each weighted by the square root of its eigenvalue.
const DISTANCE_PCS: i64 = 4;
const TOP: usize = 20;
/// Two runs on the same panel give the same axes; a PC whose present-day
/// coordinates correlate less than this between the two is a mismatch.
const MIN_AXIS_CORRELATION: f64 = 0.999;
const CAVEAT: &str = "Closeness on a PCA is genetic similarity, not descent: a genome near an ancient group \
shares ancestry with it, it does not descend from those people.";
const QPADM_NOTE: &str = "qpAdm models of ancestry proportions are not part of this step.";

/// The `.anno` columns read, matched by the start of the header cell (the
/// AADR writes long explanations after the name).
const ANNO_COLUMNS: [(Col, &str); 13] = [
    (Col::Id, "Genetic ID"),
    (Col::Individual, "Individual ID"),
    (Col::DateBp, "Date mean in BP"),
    (Col::DateSd, "Date standard deviation in BP"),
    (Col::FullDate, "Full Date"),
    (Col::Group, "Group ID"),
    (Col::Locality, "Locality"),
    (Col::Country, "Political Entity"),
    (Col::Lat, "Latitude"),
    (Col::Lon, "Longitude"),
    (Col::DataType, "Data type"),
    (
        Col::Snps,
        "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)",
    ),
    (Col::Assessment, "ASSESSMENT"),
];

/// A column of the `.anno` file this tool reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Col {
    Id,
    Individual,
    DateBp,
    DateSd,
    FullDate,
    Group,
    Locality,
    Country,
    Lat,
    Lon,
    DataType,
    Snps,
    Assessment,
}

/// One `.anno` row: the 13 cells read, stripped, in [`ANNO_COLUMNS`] order.
#[derive(Debug, Clone)]
struct Anno([String; 13]);

impl Anno {
    fn get(&self, c: Col) -> &str {
        &self.0[c as usize]
    }
}

/// A text file read line by line as Python's text mode reads it.
struct Text {
    path: PathBuf,
    lines: Lines<Box<dyn BufRead>>,
    line: String,
}

impl Text {
    fn open(path: &str, utf8: Utf8) -> Result<Text> {
        let path = PathBuf::from(path);
        let input = og_io::open(&path).map_err(|e| cannot_read(&path, &e))?;
        Ok(Text {
            path,
            lines: Lines::new(input, utf8),
            line: String::new(),
        })
    }

    /// The next line with its `\n`, or `None` at the end.
    fn raw(&mut self) -> Result<Option<&str>> {
        let more = self
            .lines
            .read_line(&mut self.line)
            .map_err(|e| cannot_read(&self.path, &e))?;
        Ok(more.then_some(self.line.as_str()))
    }

    /// The next line without its `\n`.
    fn next(&mut self) -> Result<Option<&str>> {
        Ok(self.raw()?.map(|l| l.strip_suffix('\n').unwrap_or(l)))
    }
}

/// Python's `float(s)` of a finite value, or `None`.
fn to_float(s: &str) -> Option<f64> {
    pycompat::num::parse_float(s).ok().filter(|v| v.is_finite())
}

/// Python's `float(s)`, or the input error that names the file.
fn need_float(s: &str, path: &Path) -> Result<f64> {
    pycompat::num::parse_float(s)
        .map_err(|_| Fail(format!("{}: {} is not a number", path.display(), repr(s))))
}

/// Python's `f >= i` for a float and an int: exact, whatever their sizes.
fn float_ge_int(f: f64, i: i64) -> bool {
    let fi = i as f64;
    if f != fi {
        return f > fi;
    }
    // f equals the float nearest i, so it is a whole number.
    (f as i128) >= i128::from(i)
}

/// `" ".join(str(s).replace("\t", " ").split()) or "."`.
fn clean(s: &str) -> String {
    let out = split_ws(s).collect::<Vec<_>>().join(" ");
    if out.is_empty() { ".".to_owned() } else { out }
}

// --- select -----------------------------------------------------------------------

/// The index of each read column in the `.anno` header.
fn anno_columns(header: &[&str]) -> Result<[usize; 13]> {
    let cells: Vec<&str> = header
        .iter()
        .map(|c| strip(strip(strip(c, None), Some("\"")), None))
        .collect();
    let mut out = [0usize; 13];
    for (k, (col, name)) in ANNO_COLUMNS.iter().enumerate() {
        // ASSESSMENT must not match ASSESSMENT WARNINGS: the exact cell first.
        let exact = cells.iter().position(|c| c == name);
        let hit = if *col == Col::Assessment {
            exact
        } else {
            exact.or_else(|| cells.iter().position(|c| c.starts_with(name)))
        };
        out[k] = hit.ok_or_else(|| Fail(format!("the .anno file has no column '{name}'")))?;
    }
    Ok(out)
}

fn passing(assessment: &str) -> bool {
    let a = upper(strip(assessment, None));
    PASSING.iter().any(|p| a.starts_with(p))
}

/// What [`choose`] keeps of one genetic id.
#[derive(Debug, Clone)]
enum Chosen {
    /// A present-day reference person, with their group.
    Modern(String),
    /// An ancient person to project, with their row.
    Ancient(Box<Anno>),
}

/// The present-day reference people and the ancient people to project,
/// by genetic id, one genetic id per individual (the one with the most
/// sites). Rows are streamed from the `.anno` file.
fn choose(path: &str, min_snps: i64) -> Result<HashMap<String, Chosen>> {
    let mut t = Text::open(path, Utf8::Replace)?;
    let Some(head) = t.next()? else {
        return Err(Fail(format!("{path} is empty")));
    };
    let header: Vec<&str> = head.split('\t').collect();
    let col = anno_columns(&header)?;
    let widest = *col.iter().max().expect("13 columns");
    // (kind is modern, individual) -> (sites, row), in first-seen order.
    let mut index: HashMap<(bool, String), usize> = HashMap::new();
    let mut best: Vec<(f64, Anno, bool)> = Vec::new();
    while let Some(line) = t.next()? {
        let cells: Vec<&str> = line.split('\t').collect();
        if line.is_empty() || cells.len() <= widest {
            continue;
        }
        let r = Anno(col.map(|i| strip(cells[i], None).to_owned()));
        if !passing(r.get(Col::Assessment)) {
            continue;
        }
        let bp = to_float(r.get(Col::DateBp));
        let snps = to_float(r.get(Col::Snps)).unwrap_or(0.0);
        let modern = if bp == Some(0.0)
            && MODERN_GROUPS.contains(&r.get(Col::Group))
            && r.get(Col::DataType) == MODERN_DATA_TYPE
        {
            true
        } else if bp.is_some_and(|b| b > 0.0) && float_ge_int(snps, min_snps) {
            let (Some(lat), Some(lon)) = (to_float(r.get(Col::Lat)), to_float(r.get(Col::Lon)))
            else {
                continue;
            };
            if !(LAT.0 <= lat && lat <= LAT.1 && LON.0 <= lon && lon <= LON.1) {
                continue;
            }
            false
        } else {
            continue;
        };
        let who = if r.get(Col::Individual).is_empty() {
            r.get(Col::Id)
        } else {
            r.get(Col::Individual)
        };
        let key = (modern, who.to_owned());
        match index.get(&key) {
            None => {
                index.insert(key, best.len());
                best.push((snps, r, modern));
            }
            Some(&i) if snps > best[i].0 => best[i] = (snps, r, modern),
            Some(_) => {}
        }
    }
    let mut out = HashMap::new();
    for (_, r, modern) in best {
        let id = r.get(Col::Id).to_owned();
        let c = if modern {
            Chosen::Modern(r.get(Col::Group).to_owned())
        } else {
            Chosen::Ancient(Box::new(r))
        };
        out.insert(id, c);
    }
    Ok(out)
}

const SELECT: Spec = Spec {
    verb: "select",
    usage: "usage: og ancient select --anno ANNO --ind IND --out-ind OUT_IND --poplist POPLIST --ancients ANCIENTS [--min-snps MIN_SNPS]\n",
    options: &[
        "--anno",
        "--ind",
        "--out-ind",
        "--poplist",
        "--ancients",
        "--min-snps",
    ],
    required: &["--anno", "--ind", "--out-ind", "--poplist", "--ancients"],
};

fn cmd_select(a: &Args, out: &mut dyn Write) -> Result<()> {
    let min_snps = a
        .int("--min-snps", MIN_ANCIENT_SNPS)
        .unwrap_or(MIN_ANCIENT_SNPS);
    let chosen = choose(a.req("--anno"), min_snps)?;
    let (mut n_mod, mut n_anc) = (0u64, 0u64);
    let mut groups: std::collections::BTreeSet<String> = Default::default();
    let mut ind = String::new();
    let mut t = Text::open(a.req("--ind"), Utf8::Strict)?;
    while let Some(line) = t.next()? {
        let c: Vec<&str> = split_ws(line).collect();
        if c.len() < 3 {
            continue;
        }
        let label = match chosen.get(c[0]) {
            Some(Chosen::Modern(g)) => {
                n_mod += 1;
                groups.insert(g.clone());
                g.as_str()
            }
            Some(Chosen::Ancient(_)) => {
                n_anc += 1;
                ANCIENT_LABEL
            }
            None => "Ignore",
        };
        ind += &format!("{}\t{}\t{label}\n", c[0], c[1]);
    }
    let missing: Vec<&str> = {
        let mut m: Vec<&str> = MODERN_GROUPS
            .iter()
            .copied()
            .filter(|g| !groups.contains(*g))
            .collect();
        m.sort_unstable();
        m
    };
    if n_mod < 50 || n_anc < 1 {
        return Err(Fail(format!(
            "the .ind file matches {n_mod} present-day reference people and {n_anc} ancient ones; \
             is it the .ind file of the same release as the .anno file?"
        )));
    }
    let mut files = Pending::new();
    files.put(PathBuf::from(a.req("--out-ind")), ind.as_bytes())?;
    let pops: String = groups.iter().map(|g| format!("{g}\n")).collect();
    files.put(PathBuf::from(a.req("--poplist")), pops.as_bytes())?;
    let mut anc =
        String::from("id\tgroup\tdate_bp\tdate_sd\tfull_date\tlocality\tcountry\tlat\tlon\tsnps\n");
    let sorted: BTreeMap<&String, &Chosen> = chosen.iter().collect();
    for c in sorted.values() {
        let Chosen::Ancient(r) = c else { continue };
        let cells: Vec<String> = [
            Col::Id,
            Col::Group,
            Col::DateBp,
            Col::DateSd,
            Col::FullDate,
            Col::Locality,
            Col::Country,
            Col::Lat,
            Col::Lon,
            Col::Snps,
        ]
        .iter()
        .map(|k| clean(r.get(*k)))
        .collect();
        anc += &cells.join("\t");
        anc.push('\n');
    }
    files.put(PathBuf::from(a.req("--ancients")), anc.as_bytes())?;
    files.commit_all()?;
    let suffix = if !missing.is_empty() && missing.len() <= 10 {
        format!(" (not in this file: {})", missing.join(", "))
    } else if !missing.is_empty() {
        format!(
            " ({} of the {} groups are not in this file)",
            missing.len(),
            MODERN_GROUPS.len()
        )
    } else {
        String::new()
    };
    let _ = writeln!(
        out,
        "Present-day reference people: {n_mod} in {} groups{suffix}",
        groups.len()
    );
    let _ = writeln!(
        out,
        "Ancient people projected: {n_anc} (at least {} sites, West Eurasia)",
        int_grouped(i128::from(min_snps))
    );
    Ok(())
}

// --- lift -------------------------------------------------------------------------

/// A base's complement, for the four bases.
fn complement(b: &str) -> Option<&'static str> {
    match b {
        "A" => Some("T"),
        "C" => Some("G"),
        "G" => Some("C"),
        "T" => Some("A"),
        _ => None,
    }
}

/// The chromosome column of an EIGENSOFT `.snp` line: 1 to 22, or `None`.
fn eigen_chrom(c: &str) -> Option<i64> {
    int(c).filter(|n| (1..=22).contains(n))
}

const LIFT: Spec = Spec {
    verb: "lift",
    usage: "usage: og ancient lift --snp SNP --chain CHAIN --fasta FASTA --out OUT\n",
    options: &["--snp", "--chain", "--fasta", "--out"],
    required: &["--snp", "--chain", "--fasta", "--out"],
};

fn cmd_lift(a: &Args, out: &mut dyn Write) -> Result<()> {
    let chain_path = a.req("--chain");
    let input = og_io::open(chain_path).map_err(|e| cannot_read(Path::new(chain_path), &e))?;
    let chain = Chain::read(input).map_err(|e| Fail(format!("{chain_path}: {e}")))?;
    let fasta_path = a.req("--fasta");
    let mut fasta = Indexed::open(fasta_path).map_err(|e| Fail(format!("{fasta_path}: {e}")))?;
    let mut kept: Vec<(i64, i64, String, &'static str, &'static str)> = Vec::new();
    let mut dropped: BTreeMap<&str, u64> = BTreeMap::new();
    let mut n: u64 = 0;
    let snp_path = a.req("--snp");
    let mut t = Text::open(snp_path, Utf8::Replace)?;
    let mut lineno = 0u64;
    while let Some(line) = t.next()? {
        lineno += 1;
        let c: Vec<&str> = split_ws(line).collect();
        if c.len() < 6 {
            continue;
        }
        n += 1;
        let Some(chrom) = eigen_chrom(c[1]) else {
            *dropped.entry("not an autosome").or_default() += 1;
            continue;
        };
        let (u1, u2) = (upper(c[4]), upper(c[5]));
        let (Some(k1), Some(k2)) = (complement(&u1), complement(&u2)) else {
            *dropped
                .entry("alleles not two different bases")
                .or_default() += 1;
            continue;
        };
        if u1 == u2 {
            *dropped
                .entry("alleles not two different bases")
                .or_default() += 1;
            continue;
        }
        let pos = pycompat::num::parse_float(c[3])
            .ok()
            .filter(|p| p.is_finite())
            .map(f64::trunc)
            .filter(|p| p.abs() < 9.0e18)
            .ok_or_else(|| {
                Fail(format!(
                    "{snp_path}: line {lineno}: position {} is not a number",
                    repr(c[3])
                ))
            })? as i64;
        let name = format!("chr{chrom}");
        let Some(hit) = chain.lift(&name, pos) else {
            *dropped.entry("no single place on GRCh38").or_default() += 1;
            continue;
        };
        if hit.chrom != name {
            *dropped.entry("moves to another chromosome").or_default() += 1;
            continue;
        }
        let (a1, a2) = if hit.strand == "-" {
            (k1, k2)
        } else {
            (
                complement(k1).expect("a base"),
                complement(k2).expect("a base"),
            )
        };
        let q_pos = hit.pos;
        let base = if q_pos >= 1 {
            fasta
                .base(hit.chrom, (q_pos - 1) as u64)
                .map_err(|e| Fail(format!("{fasta_path}: {e}")))?
        } else {
            None
        };
        let base = base.map(|b| b.to_ascii_uppercase());
        if !base.is_some_and(|b| b == a1.as_bytes()[0] || b == a2.as_bytes()[0]) {
            *dropped
                .entry("neither allele is the GRCh38 base")
                .or_default() += 1;
            continue;
        }
        kept.push((chrom, q_pos, c[0].to_owned(), a1, a2));
    }
    kept.sort_unstable();
    let mut text = String::from("snp\tchrom\tpos\tallele1\tallele2\n");
    for (chrom, pos, sid, a1, a2) in &kept {
        text += &format!("{sid}\tchr{chrom}\t{pos}\t{a1}\t{a2}\n");
    }
    let mut files = Pending::new();
    files.put(PathBuf::from(a.req("--out")), text.as_bytes())?;
    files.commit_all()?;
    let _ = writeln!(
        out,
        "Sites of the .snp file: {}; on GRCh38: {}",
        int_grouped(i128::from(n)),
        int_grouped(kept.len() as i128)
    );
    for (why, k) in &dropped {
        let _ = writeln!(out, "  dropped, {why}: {}", int_grouped(i128::from(*k)));
    }
    if kept.is_empty() {
        return Err(Fail(
            "no site could be placed on GRCh38: check the chain file and the reference".into(),
        ));
    }
    Ok(())
}

// --- person-geno ------------------------------------------------------------------

/// Copies of allele1 in a bcftools `%TGT` genotype, or `None` (missing, a
/// haploid or half call, or an allele that is neither of the two).
fn dosage(gt: &str, a1: &str, a2: &str) -> Option<u8> {
    let mut parts = gt.split(['/', '|']);
    let (Some(x), Some(y), None) = (parts.next(), parts.next(), parts.next()) else {
        return None;
    };
    let mut d = 0;
    for al in [x, y] {
        if al == a1 {
            d += 1;
        } else if al != a2 {
            return None;
        }
    }
    Some(d)
}

/// The panel's GRCh38 sites: per (chromosome, position), the sites there
/// (id, allele1, allele2), packed in one array.
struct Sites {
    chroms: HashMap<String, u32>,
    at: HashMap<(u32, i64), Vec<u32>>,
    sites: Vec<(Box<str>, Box<str>, Box<str>)>,
}

fn read_sites(path: &str) -> Result<Sites> {
    let mut s = Sites {
        chroms: HashMap::new(),
        at: HashMap::new(),
        sites: Vec::new(),
    };
    let mut t = Text::open(path, Utf8::Strict)?;
    let _ = t.next()?;
    let mut lineno = 1u64;
    while let Some(line) = t.next()? {
        lineno += 1;
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            continue;
        }
        let pos = int(c[2]).ok_or_else(|| {
            Fail(format!(
                "{path}: line {lineno}: position {} is not a whole number",
                repr(c[2])
            ))
        })?;
        let next = s.chroms.len() as u32;
        let chrom = *s.chroms.entry(c[1].to_owned()).or_insert(next);
        s.at.entry((chrom, pos))
            .or_default()
            .push(s.sites.len() as u32);
        s.sites.push((c[0].into(), c[3].into(), c[4].into()));
    }
    Ok(s)
}

const PERSON_GENO: Spec = Spec {
    verb: "person-geno",
    usage: "usage: og ancient person-geno --sites SITES --genotypes GENOTYPES --snp SNP --moderns-geno MODERNS_GENO --out OUT --stats STATS\n",
    options: &[
        "--sites",
        "--genotypes",
        "--snp",
        "--moderns-geno",
        "--out",
        "--stats",
    ],
    required: &[
        "--sites",
        "--genotypes",
        "--snp",
        "--moderns-geno",
        "--out",
        "--stats",
    ],
};

fn cmd_person_geno(a: &Args, out: &mut dyn Write) -> Result<()> {
    let sites = read_sites(a.req("--sites"))?;
    // The person's dosage by site id: a later line for the same site wins.
    let mut dos: HashMap<&str, u8> = HashMap::new();
    let (mut called, mut other) = (0u64, 0u64);
    let gpath = a.req("--genotypes");
    let mut t = Text::open(gpath, Utf8::Strict)?;
    let mut lineno = 0u64;
    while let Some(line) = t.next()? {
        lineno += 1;
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            continue;
        }
        let pos = int(c[1]).ok_or_else(|| {
            Fail(format!(
                "{gpath}: line {lineno}: position {} is not a whole number",
                repr(c[1])
            ))
        })?;
        let Some(&chrom) = sites.chroms.get(c[0]) else {
            continue;
        };
        let Some(here) = sites.at.get(&(chrom, pos)) else {
            continue;
        };
        let gt = upper(c[4]);
        // quirk: a position the AADR holds two SNPs at is counted once per
        // SNP in sites_genotyped, as the Python counted it.
        for &i in here {
            let (sid, a1, a2) = &sites.sites[i as usize];
            match dosage(&gt, a1, a2) {
                None => other += 1,
                Some(d) => {
                    dos.insert(sid, d);
                    called += 1;
                }
            }
        }
    }
    let snp_path = a.req("--snp");
    let geno_path = a.req("--moderns-geno");
    let mut snp = Text::open(snp_path, Utf8::Strict)?;
    let mut geno = Text::open(geno_path, Utf8::Strict)?;
    let mut files = Pending::new();
    let out_path = PathBuf::from(a.req("--out"));
    let shown = out_path.clone();
    let o = files.create(out_path)?;
    let (mut n, mut used) = (0u64, 0u64);
    let mut width: Option<usize> = None;
    let mut row = String::new();
    while let Some(sline) = snp.next()? {
        let Some(id) = split_ws(sline).next() else {
            continue;
        };
        let d = dos.get(id).copied();
        let g = geno.next()?.unwrap_or("");
        if g.is_empty() {
            return Err(Fail(format!("{geno_path} has fewer lines than {snp_path}")));
        }
        let w = g.chars().count();
        match width {
            None => width = Some(w),
            Some(want) if w != want => {
                return Err(Fail(format!(
                    "{geno_path} line {} has {w} genotypes, not {want}",
                    n + 1
                )));
            }
            Some(_) => {}
        }
        if d.is_some() {
            used += 1;
        }
        row.clear();
        row.push_str(g);
        match d {
            None => row.push('9'),
            Some(d) => row.push(char::from(b'0' + d)),
        }
        row.push('\n');
        o.write_all(row.as_bytes())
            .map_err(|e| super::cannot_write(&shown, &e))?;
        n += 1;
    }
    if geno.raw()?.is_some() {
        return Err(Fail(format!("{geno_path} has more lines than {snp_path}")));
    }
    let genotyped = called + other;
    files.put(
        PathBuf::from(a.req("--stats")),
        format!(
            "panel_sites\t{n}\nsites_genotyped\t{genotyped}\nother_allele\t{other}\nsites_used\t{used}\n"
        )
        .as_bytes(),
    )?;
    files.commit_all()?;
    let _ = writeln!(
        out,
        "Sites of the panel: {}; genotyped from the gVCF: {}; used: {} (an allele that is neither of a site's two: {})",
        int_grouped(i128::from(n)),
        int_grouped(i128::from(genotyped)),
        int_grouped(i128::from(used)),
        int_grouped(i128::from(other))
    );
    if used == 0 {
        return Err(Fail(
            "none of the panel's sites was genotyped: is the gVCF on GRCh38 with chr names?".into(),
        ));
    }
    Ok(())
}

// --- evec, align ------------------------------------------------------------------

/// A smartpca `.evec`: the eigenvalues, and each row's label and PCs in
/// file order (a repeated id keeps its first place and its last values).
struct Evec {
    evals: Vec<f64>,
    rows: Vec<(String, String, Vec<f64>)>,
    index: HashMap<String, usize>,
}

impl Evec {
    fn get(&self, id: &str) -> Option<(&str, &[f64])> {
        self.index.get(id).map(|&i| {
            let r = &self.rows[i];
            (r.1.as_str(), r.2.as_slice())
        })
    }
}

fn read_evec(path: &str) -> Result<Evec> {
    let mut e = Evec {
        evals: Vec::new(),
        rows: Vec::new(),
        index: HashMap::new(),
    };
    let p = Path::new(path);
    let mut t = Text::open(path, Utf8::Strict)?;
    while let Some(line) = t.next()? {
        let c: Vec<&str> = split_ws(line).collect();
        let Some(first) = c.first() else { continue };
        if first.starts_with('#') {
            if *first == "#eigvals:" {
                e.evals = c[1..]
                    .iter()
                    .map(|x| need_float(x, p))
                    .collect::<Result<_>>()?;
            }
            continue;
        }
        let label = c[c.len() - 1].to_owned();
        let pcs: Vec<f64> = if c.len() > 1 {
            c[1..c.len() - 1]
                .iter()
                .map(|x| need_float(x, p))
                .collect::<Result<_>>()?
        } else {
            Vec::new()
        };
        match e.index.get(*first) {
            Some(&i) => e.rows[i] = ((*first).to_owned(), label, pcs),
            None => {
                e.index.insert((*first).to_owned(), e.rows.len());
                e.rows.push(((*first).to_owned(), label, pcs));
            }
        }
    }
    if e.evals.is_empty() || e.rows.is_empty() {
        return Err(Fail(format!(
            "{path} is not a smartpca .evec file (no #eigvals: line or no rows)"
        )));
    }
    Ok(e)
}

/// The Pearson correlation of two equal-length series, 0 when either is
/// flat; means by `fmean`, sums by the builtin `sum`.
fn correlation(x: &[f64], y: &[f64]) -> f64 {
    let mx = pycompat::stats::fmean(x.iter().copied()).unwrap_or(f64::NAN);
    let my = pycompat::stats::fmean(y.iter().copied()).unwrap_or(f64::NAN);
    let sxy = pycompat::num::sum(x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)));
    let sx = pycompat::num::sum(x.iter().map(|a| (a - mx) * (a - mx))).sqrt();
    let sy = pycompat::num::sum(y.iter().map(|b| (b - my) * (b - my))).sqrt();
    if sx > 0.0 && sy > 0.0 {
        sxy / (sx * sy)
    } else {
        0.0
    }
}

/// The k-th PC of a row, or the error that names the row.
fn pc(pcs: &[f64], k: usize, id: &str, path: &str) -> Result<f64> {
    pcs.get(k)
        .copied()
        .ok_or_else(|| Fail(format!("{path}: {id} has no PC{}", k + 1)))
}

/// The sign that turns each of the run's PCs into the panel's, from the
/// present-day people both runs used to build the axes.
fn axis_signs(
    panel: &Evec,
    panel_path: &str,
    run: &Evec,
    run_path: &str,
    n_pcs: usize,
) -> Result<Vec<f64>> {
    let mut shared: Vec<&str> = run
        .rows
        .iter()
        .filter(|(id, lab, _)| {
            lab != ANCIENT_LABEL
                && lab != PERSON_LABEL
                && panel.get(id).is_some_and(|(pl, _)| pl == lab)
        })
        .map(|(id, ..)| id.as_str())
        .collect();
    shared.sort_unstable();
    if shared.len() < 10 {
        return Err(Fail(format!(
            "only {} present-day people are in both runs; the panel and this run do not match",
            shared.len()
        )));
    }
    let mut signs = Vec::with_capacity(n_pcs);
    for k in 0..n_pcs {
        let mut x = Vec::with_capacity(shared.len());
        let mut y = Vec::with_capacity(shared.len());
        for id in &shared {
            x.push(pc(panel.get(id).expect("shared").1, k, id, panel_path)?);
            y.push(pc(run.get(id).expect("shared").1, k, id, run_path)?);
        }
        let r = correlation(&x, &y);
        if r.abs() < MIN_AXIS_CORRELATION {
            return Err(Fail(format!(
                "PC{} of this run correlates {} with the panel's over {} present-day people: \
                 the axes differ. Install the panel again (setup.sh --ancient-panel).",
                k + 1,
                fixed(r, 4),
                shared.len()
            )));
        }
        signs.push(if r > 0.0 { 1.0 } else { -1.0 });
    }
    Ok(signs)
}

const ALIGN: Spec = Spec {
    verb: "align",
    usage: "usage: og ancient align --panel-evec PANEL_EVEC --run-evec RUN_EVEC [--pcs PCS]\n",
    options: &["--panel-evec", "--run-evec", "--pcs"],
    required: &["--panel-evec", "--run-evec"],
};

fn cmd_align(a: &Args, out: &mut dyn Write) -> Result<()> {
    let (pp, rp) = (a.req("--panel-evec"), a.req("--run-evec"));
    let panel = read_evec(pp)?;
    let run = read_evec(rp)?;
    let pcs = a.int("--pcs", DISTANCE_PCS).unwrap_or(DISTANCE_PCS);
    let signs = axis_signs(&panel, pp, &run, rp, usize::try_from(pcs).unwrap_or(0))?;
    let words: Vec<String> = signs
        .iter()
        .enumerate()
        .map(|(i, s)| format!("PC{}{}", i + 1, if *s > 0.0 { '+' } else { '-' }))
        .collect();
    let _ = writeln!(out, "Axes agree with the panel's: {}", words.join(" "));
    Ok(())
}

// --- report -----------------------------------------------------------------------

/// A whole float as Python prints the int of it.
fn whole(x: f64) -> String {
    int_str(x)
}

/// `int(round(x))` printed.
fn rounded(x: f64) -> String {
    whole(round_int(x).unwrap_or(f64::NAN))
}

/// A date in years before 1950 as `N BCE` or `N CE`, the AADR's own
/// reading: BCE = BP - 1950 and CE = 1950 - BP, with 1950 BP as 1 BCE.
fn year_label(bp: f64) -> String {
    let y = 1950.0 - round_int(bp).unwrap_or(0.0);
    if y >= 1.0 {
        format!("{} CE", whole(y))
    } else {
        format!("{} BCE", whole(if -y > 1.0 { -y } else { 1.0 }))
    }
}

/// The span of a population's dates, oldest first: `5300-4900 BCE`,
/// `200-450 CE`, or `300 BCE-200 CE`.
fn range_label(old_bp: f64, young_bp: f64) -> String {
    let (a, b) = (year_label(old_bp), year_label(young_bp));
    if a == b {
        return a;
    }
    let (na, ea) = a.split_once(' ').expect("N ERA");
    let (nb, eb) = b.split_once(' ').expect("N ERA");
    if ea == eb {
        format!("{na}-{nb} {ea}")
    } else {
        format!("{a}-{b}")
    }
}

/// An `ancients.tsv` row as `csv.DictReader` gives it: a short row's
/// missing cells are `None`.
struct AncientRow {
    group: String,
    date_bp: Option<String>,
    country: Option<String>,
    locality: Option<String>,
}

/// `ancients.tsv` by id, in file order (a repeated id keeps its first
/// place and its last values).
fn read_ancients(path: &str) -> Result<Vec<(String, AncientRow)>> {
    let mut t = Text::open(path, Utf8::Strict)?;
    let mut rows: Vec<(String, AncientRow)> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    // DictReader: the first non-empty row is the header.
    let mut names: Option<Vec<String>> = None;
    while let Some(line) = t.next()? {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let Some(names) = names.as_ref() else {
            names = Some(cells.iter().map(|s| (*s).to_owned()).collect());
            continue;
        };
        // The last column of a repeated name wins, as in a dict.
        let field = |key: &str| -> Option<String> {
            let k = names.iter().rposition(|n| n == key)?;
            cells.get(k).map(|s| (*s).to_owned())
        };
        if !names.iter().any(|n| n == "id") {
            return Err(Fail(format!("{path} has no id column")));
        }
        let id = field("id").unwrap_or_else(|| "None".to_owned());
        let r = AncientRow {
            group: field("group").unwrap_or_else(|| "None".to_owned()),
            date_bp: field("date_bp"),
            country: field("country"),
            locality: field("locality"),
        };
        match index.get(&id) {
            Some(&i) => rows[i].1 = r,
            None => {
                index.insert(id.clone(), rows.len());
                rows.push((id, r));
            }
        }
    }
    Ok(rows)
}

/// The most common non-empty value other than `.`, the first in sort
/// order on a tie; `.` when there is none.
fn most_common<'a>(values: impl Iterator<Item = Option<&'a str>>) -> String {
    let mut counts: BTreeMap<&str, u64> = BTreeMap::new();
    for v in values.flatten() {
        if !v.is_empty() && v != "." {
            *counts.entry(v).or_default() += 1;
        }
    }
    let top = counts.values().copied().max().unwrap_or(0);
    counts
        .into_iter()
        .find(|(_, k)| *k == top)
        .map_or_else(|| ".".to_owned(), |(v, _)| v.to_owned())
}

/// One ancient group, nearest first.
struct Near {
    rank: usize,
    group: String,
    n: usize,
    distance: f64,
    pc1: f64,
    pc2: f64,
    date_bp: Option<f64>,
    dates: String,
    country: String,
    locality: String,
}

/// The ancient groups ranked by the weighted distance from the person to
/// the group's mean position on the first `n_pcs` PCs; the top [`TOP`] and
/// the number of groups.
fn nearest_populations(
    person: &[f64],
    panel: &Evec,
    ancients: &[(String, AncientRow)],
    n_pcs: usize,
) -> Result<(Vec<Near>, usize)> {
    let w: Vec<f64> = panel
        .evals
        .iter()
        .take(n_pcs)
        .map(|v| if 0.0 > *v { 0.0 } else { *v }.sqrt())
        .collect();
    // Groups in first-seen order, each with its members.
    let mut order: Vec<&str> = Vec::new();
    let mut groups: HashMap<&str, Vec<(&[f64], &AncientRow)>> = HashMap::new();
    for (gid, r) in ancients {
        if let Some((label, pcs)) = panel.get(gid)
            && label == ANCIENT_LABEL
        {
            let g = groups.entry(r.group.as_str()).or_insert_with(|| {
                order.push(r.group.as_str());
                Vec::new()
            });
            g.push((pcs, r));
        }
    }
    let mut rows = Vec::with_capacity(order.len());
    for g in order {
        let members = &groups[g];
        let k = members[0].0.len().min(n_pcs);
        let mut centre = Vec::with_capacity(k);
        for i in 0..k {
            let mut xs = Vec::with_capacity(members.len());
            for (pcs, _) in members {
                xs.push(*pcs.get(i).ok_or_else(|| {
                    Fail(format!(
                        "a member of the ancient group {g} has no PC{}",
                        i + 1
                    ))
                })?);
            }
            centre.push(pycompat::stats::fmean(xs).map_err(|e| Fail(e.to_string()))?);
        }
        let Some(&pc1) = centre.first() else {
            return Err(Fail(format!("the ancient group {g} has no PCs")));
        };
        let distance = pycompat::num::sum(person.iter().take(n_pcs).zip(&centre).zip(&w).map(
            |((a, b), wi)| {
                let d = wi * (a - b);
                d * d
            },
        ))
        .sqrt();
        let dates: Vec<f64> = members
            .iter()
            .filter_map(|(_, r)| r.date_bp.as_deref().and_then(to_float))
            .collect();
        let (date_bp, dates_text) = if dates.is_empty() {
            (None, ".".to_owned())
        } else {
            let old = dates
                .iter()
                .copied()
                .fold(dates[0], |m, d| if d > m { d } else { m });
            let young = dates
                .iter()
                .copied()
                .fold(dates[0], |m, d| if d < m { d } else { m });
            (
                Some(
                    pycompat::stats::median(dates.iter().copied())
                        .map_err(|e| Fail(e.to_string()))?,
                ),
                range_label(old, young),
            )
        };
        rows.push(Near {
            rank: 0,
            group: g.to_owned(),
            n: members.len(),
            distance,
            pc1,
            pc2: centre.get(1).copied().unwrap_or(0.0),
            date_bp,
            dates: dates_text,
            country: most_common(members.iter().map(|(_, r)| r.country.as_deref())),
            locality: most_common(members.iter().map(|(_, r)| r.locality.as_deref())),
        });
    }
    rows.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.group.cmp(&b.group))
    });
    for (i, r) in rows.iter_mut().enumerate() {
        r.rank = i + 1;
    }
    let n = rows.len();
    rows.truncate(TOP);
    Ok((rows, n))
}

/// Age bins of the plot: (oldest BP, colour, legend), oldest first.
const AGE_BINS: [(f64, &str, &str); 4] = [
    (f64::INFINITY, "#5aa9e6", "before 8000 BCE"),
    (9950.0, "#7fc97f", "8000-4000 BCE"),
    (5950.0, "#f0a35e", "4000-1000 BCE"),
    (2950.0, "#e377c2", "after 1000 BCE"),
];

fn age_colour(bp: Option<f64>) -> &'static str {
    let mut colour = AGE_BINS[0].1;
    for (limit, c, _) in AGE_BINS {
        if bp.is_some_and(|b| b < limit) {
            colour = c;
        }
    }
    colour
}

/// Python's `min()` and `max()` of a non-empty list of floats.
fn py_min_max(xs: &[f64]) -> (f64, f64) {
    let (mut lo, mut hi) = (xs[0], xs[0]);
    for &x in &xs[1..] {
        if x < lo {
            lo = x;
        }
        if x > hi {
            hi = x;
        }
    }
    (lo, hi)
}

/// The PC1/PC2 plane: present-day people grey, ancient people coloured by
/// age, the person a large ring, on a dark background.
fn svg_plot(
    person: &[f64],
    panel: &Evec,
    ancients: &HashMap<String, &AncientRow>,
    sample: &str,
    nearest: &[Near],
) -> Result<String> {
    const WIDTH: i64 = 900;
    const HEIGHT: i64 = 720;
    const PAD: f64 = 0.04;
    const LEFT: i64 = 60;
    const RIGHT: i64 = 20;
    const TOPM: i64 = 50;
    const BOTTOM: i64 = 110;
    let two = |id: &str, pcs: &[f64]| -> Result<(f64, f64)> {
        match pcs {
            [x, y, ..] => Ok((*x, *y)),
            _ => Err(Fail(format!("{id} has fewer than two PCs to plot"))),
        }
    };
    let mut xs = Vec::with_capacity(panel.rows.len() + 1);
    let mut ys = Vec::with_capacity(panel.rows.len() + 1);
    for (id, _, pcs) in &panel.rows {
        let (x, y) = two(id, pcs)?;
        xs.push(x);
        ys.push(y);
    }
    let (px0, py0) = two(sample, person)?;
    xs.push(px0);
    ys.push(py0);
    let (x0, x1) = py_min_max(&xs);
    let (y0, y1) = py_min_max(&ys);
    let or_one = |d: f64| if d == 0.0 { 1.0 } else { d };
    let (dx, dy) = (or_one(x1 - x0), or_one(y1 - y0));
    let (pw, ph) = (WIDTH - LEFT - RIGHT, HEIGHT - TOPM - BOTTOM);
    let (left, top, pwf, phf) = (LEFT as f64, TOPM as f64, pw as f64, ph as f64);
    let x = |v: f64| left + ((v - x0) / dx * (1.0 - 2.0 * PAD) + PAD) * pwf;
    let y = |v: f64| top + (1.0 - ((v - y0) / dy * (1.0 - 2.0 * PAD) + PAD)) * phf;
    let e = |s: &str| html_escape(s).into_owned();
    let mut out = vec![
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{WIDTH}\" height=\"{HEIGHT}\" viewBox=\"0 0 {WIDTH} {HEIGHT}\" font-family=\"sans-serif\" font-size=\"13\">"
        ),
        "<style>.m{fill:#9a9a9a;fill-opacity:.55}.a{fill-opacity:.7}text{fill:#e6e6e6}</style>"
            .to_owned(),
        format!("<rect width=\"{WIDTH}\" height=\"{HEIGHT}\" fill=\"#15171c\"/>"),
        format!(
            "<text x=\"{LEFT}\" y=\"28\" font-size=\"16\">{} among ancient West Eurasians (AADR, projected onto present-day PCs)</text>",
            e(sample)
        ),
        format!(
            "<rect x=\"{LEFT}\" y=\"{TOPM}\" width=\"{pw}\" height=\"{ph}\" fill=\"none\" stroke=\"#444\"/>"
        ),
    ];
    for (_, label, pcs) in &panel.rows {
        if label != ANCIENT_LABEL {
            out.push(format!(
                "<circle class=\"m\" cx=\"{}\" cy=\"{}\" r=\"1.6\"/>",
                fixed(x(pcs[0]), 1),
                fixed(y(pcs[1]), 1)
            ));
        }
    }
    for (gid, label, pcs) in &panel.rows {
        if label == ANCIENT_LABEL
            && let Some(r) = ancients.get(gid)
        {
            let c = age_colour(r.date_bp.as_deref().and_then(to_float));
            out.push(format!(
                "<circle class=\"a\" cx=\"{}\" cy=\"{}\" r=\"2\" fill=\"{c}\"/>",
                fixed(x(pcs[0]), 1),
                fixed(y(pcs[1]), 1)
            ));
        }
    }
    for r in nearest.iter().take(5) {
        out.push(format!(
            "<text x=\"{}\" y=\"{}\" font-size=\"11\">{}</text>",
            fixed(x(r.pc1) + 6.0, 1),
            fixed(y(r.pc2) - 4.0, 1),
            e(&r.group)
        ));
    }
    let (px, py) = (x(px0), y(py0));
    out.push(format!(
        "<circle cx=\"{}\" cy=\"{}\" r=\"9\" fill=\"none\" stroke=\"#ffffff\" stroke-width=\"3\"/>",
        fixed(px, 1),
        fixed(py, 1)
    ));
    out.push(format!(
        "<circle cx=\"{}\" cy=\"{}\" r=\"3\" fill=\"#ffd400\"/>",
        fixed(px, 1),
        fixed(py, 1)
    ));
    out.push(format!(
        "<text x=\"{}\" y=\"{}\" font-weight=\"bold\">{}</text>",
        fixed(px + 12.0, 1),
        fixed(py + 4.0, 1),
        e(sample)
    ));
    out.push(format!(
        "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\">PC1</text>",
        fixed(left + pwf / 2.0, 0),
        TOPM + ph + 22
    ));
    let mid = fixed(top + phf / 2.0, 0);
    out.push(format!(
        "<text x=\"18\" y=\"{mid}\" transform=\"rotate(-90 18 {mid})\" text-anchor=\"middle\">PC2</text>"
    ));
    let (mut lx, ly) = (LEFT, TOPM + ph + 48);
    out.push(format!(
        "<circle class=\"m\" cx=\"{}\" cy=\"{}\" r=\"4\"/><text x=\"{}\" y=\"{ly}\">present day</text>",
        lx + 5,
        ly - 4,
        lx + 14
    ));
    lx += 120;
    for (_, c, name) in AGE_BINS {
        out.push(format!(
            "<circle cx=\"{}\" cy=\"{}\" r=\"4\" fill=\"{c}\"/><text x=\"{}\" y=\"{ly}\">{}</text>",
            lx + 5,
            ly - 4,
            lx + 14,
            e(name)
        ));
        lx += 150;
    }
    out.push(format!(
        "<text x=\"{LEFT}\" y=\"{}\" font-size=\"12\">{}</text>",
        ly + 26,
        e(CAVEAT)
    ));
    out.push("</svg>".to_owned());
    Ok(out.join("\n") + "\n")
}

const REPORT: Spec = Spec {
    verb: "report",
    usage: "usage: og ancient report --panel-evec PANEL_EVEC --run-evec RUN_EVEC --ancients ANCIENTS --sample SAMPLE --release RELEASE [--stats STATS] --out-tsv OUT_TSV --out-info OUT_INFO --out-txt OUT_TXT --out-svg OUT_SVG\n",
    options: &[
        "--panel-evec",
        "--run-evec",
        "--ancients",
        "--sample",
        "--release",
        "--stats",
        "--out-tsv",
        "--out-info",
        "--out-txt",
        "--out-svg",
    ],
    required: &[
        "--panel-evec",
        "--run-evec",
        "--ancients",
        "--sample",
        "--release",
        "--out-tsv",
        "--out-info",
        "--out-txt",
        "--out-svg",
    ],
};

fn cmd_report(a: &Args, out: &mut dyn Write) -> Result<()> {
    let (pp, rp) = (a.req("--panel-evec"), a.req("--run-evec"));
    let panel = read_evec(pp)?;
    let run = read_evec(rp)?;
    let sample = a.req("--sample");
    let Some((PERSON_LABEL, pcs)) = run.get(sample) else {
        return Err(Fail(format!(
            "{rp} has no row for {sample} labelled {PERSON_LABEL}"
        )));
    };
    let n_pcs = (DISTANCE_PCS as usize)
        .min(panel.evals.len())
        .min(pcs.len());
    let signs = axis_signs(&panel, pp, &run, rp, n_pcs.max(2))?;
    let person: Vec<f64> = signs.iter().zip(pcs).map(|(s, v)| s * v).collect();
    let ancients = read_ancients(a.req("--ancients"))?;
    let (nearest, n_groups) = nearest_populations(&person, &panel, &ancients, n_pcs)?;
    let mut stats: HashMap<String, String> = HashMap::new();
    if let Some(sp) = a.get("--stats").filter(|s| !s.is_empty()) {
        let mut t = Text::open(sp, Utf8::Strict)?;
        while let Some(line) = t.next()? {
            let (k, v) = line.split_once('\t').unwrap_or((line, ""));
            stats.insert(k.to_owned(), v.to_owned());
        }
    }
    let stat = |k: &str| stats.get(k).map_or(".", String::as_str).to_owned();
    let n_anc = panel.rows.iter().filter(|r| r.1 == ANCIENT_LABEL).count();
    let n_mod = panel
        .rows
        .iter()
        .filter(|r| r.1 != ANCIENT_LABEL && r.1 != PERSON_LABEL)
        .count();
    let mut files = Pending::new();
    let mut tsv =
        String::from("rank\tgroup\tn\tdistance\tpc1\tpc2\tdate_bp\tdates\tcountry\tlocality\n");
    for r in &nearest {
        tsv += &format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.rank,
            clean(&r.group),
            r.n,
            fixed(r.distance, 5),
            fixed(r.pc1, 5),
            fixed(r.pc2, 5),
            r.date_bp.map_or_else(|| ".".to_owned(), rounded),
            r.dates,
            clean(&r.country),
            clean(&r.locality)
        );
    }
    files.put(PathBuf::from(a.req("--out-tsv")), tsv.as_bytes())?;
    let release = a.req("--release");
    let p1 = person.first().copied().unwrap_or(f64::NAN);
    let p2 = person.get(1).copied().unwrap_or(f64::NAN);
    let info: [(&str, String); 9] = [
        ("release", release.to_owned()),
        ("sites_used", stat("sites_used")),
        ("panel_sites", stat("panel_sites")),
        ("present_day_people", n_mod.to_string()),
        ("ancient_people", n_anc.to_string()),
        ("ancient_groups", n_groups.to_string()),
        ("pcs_in_distance", n_pcs.to_string()),
        ("pc1", fixed(p1, 5)),
        ("pc2", fixed(p2, 5)),
    ];
    let mut text = String::from("key\tvalue\n");
    for (k, v) in &info {
        text += &format!("{k}\t{v}\n");
    }
    files.put(PathBuf::from(a.req("--out-info")), text.as_bytes())?;
    let by_id: HashMap<String, &AncientRow> =
        ancients.iter().map(|(id, r)| (id.clone(), r)).collect();
    let svg = svg_plot(&person, &panel, &by_id, sample, &nearest)?;
    let mut lines = vec![
        format!(
            "Ancient populations nearest to {sample} (AADR {release}, smartpca with lsqproject)"
        ),
        format!(
            "Sites used: {} of {}; axes from {n_mod} present-day West Eurasians; {n_anc} ancient people in {n_groups} groups projected.",
            stat("sites_used"),
            stat("panel_sites")
        ),
        format!(
            "Distance: the first {n_pcs} PCs, each weighted by the square root of its eigenvalue."
        ),
        String::new(),
    ];
    for r in nearest.iter().take(10) {
        lines.push(format!(
            "  {}. {} {} {} n={} d={}",
            rjust(&r.rank.to_string(), 2),
            ljust(&r.group, 45),
            ljust(&r.dates, 22),
            ljust(&r.country, 20),
            ljust(&r.n.to_string(), 4),
            fixed(r.distance, 4)
        ));
    }
    lines.extend([String::new(), CAVEAT.to_owned(), QPADM_NOTE.to_owned()]);
    let body = lines.join("\n") + "\n";
    files.put(PathBuf::from(a.req("--out-svg")), svg.as_bytes())?;
    files.put(PathBuf::from(a.req("--out-txt")), body.as_bytes())?;
    files.commit_all()?;
    let _ = out.write_all(body.as_bytes());
    Ok(())
}

// --- the command line -------------------------------------------------------------

const USAGE: &str = "usage: og ancient {select,lift,person-geno,align,report} ...\n";

type Run = fn(&Args, &mut dyn Write) -> Result<()>;

/// `og ancient <verb> [args]`: the exit code. Results go to `out`,
/// messages to `err`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let Some(verb) = args.first() else {
        let _ = write!(err, "{USAGE}");
        let _ = writeln!(
            err,
            "og ancient: error: the following arguments are required: verb"
        );
        return 2;
    };
    let (spec, run): (&Spec, Run) = match verb.as_str() {
        "select" => (&SELECT, cmd_select),
        "lift" => (&LIFT, cmd_lift),
        "person-geno" => (&PERSON_GENO, cmd_person_geno),
        "align" => (&ALIGN, cmd_align),
        "report" => (&REPORT, cmd_report),
        "-h" | "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        other => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(
                err,
                "og ancient: error: argument verb: invalid choice: {}",
                repr(other)
            );
            return 2;
        }
    };
    let parsed = Args::parse(spec, &args[1..]).and_then(|a| {
        for name in ["--min-snps", "--pcs"] {
            if spec.options.contains(&name) {
                a.int(name, 0)?;
            }
        }
        Ok(a)
    });
    let a = match parsed {
        Ok(a) => a,
        Err(ArgError::Help) => {
            let _ = write!(out, "{}", spec.usage);
            return 0;
        }
        Err(ArgError::Usage(msg)) => {
            let _ = write!(err, "{}", spec.usage);
            let _ = writeln!(err, "og ancient {}: error: {msg}", spec.verb);
            return 2;
        }
    };
    match run(&a, out) {
        Ok(()) => 0,
        Err(Fail(msg)) => {
            let _ = out.flush();
            let _ = writeln!(err, "ERROR: {msg}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dosages() {
        assert_eq!(dosage("A/G", "A", "G"), Some(1));
        assert_eq!(dosage("A|A", "A", "G"), Some(2));
        assert_eq!(dosage("G/G", "A", "G"), Some(0));
        assert_eq!(dosage("A/T", "A", "G"), None);
        assert_eq!(dosage("./.", "A", "G"), None);
        assert_eq!(dosage("A", "A", "G"), None);
    }

    #[test]
    fn year_labels() {
        assert_eq!(year_label(4250.0), "2300 BCE");
        assert_eq!(year_label(1950.0), "1 BCE");
        assert_eq!(year_label(1949.0), "1 CE");
        assert_eq!(year_label(1750.5), "200 CE");
        assert_eq!(year_label(8025.0), "6075 BCE");
        assert_eq!(year_label(450.0), "1500 CE");
        assert_eq!(range_label(3000.0, 3000.0), "1050 BCE");
        assert_eq!(range_label(7250.0, 6850.0), "5300-4900 BCE");
        assert_eq!(range_label(1750.0, 1500.0), "200-450 CE");
        assert_eq!(range_label(2250.0, 1750.0), "300 BCE-200 CE");
        assert_eq!(range_label(2250.0, 2250.0), "300 BCE");
    }

    #[test]
    fn most_common_breaks_ties_by_name() {
        let v = [
            Some("b"),
            Some("a"),
            Some("b"),
            Some("a"),
            Some("."),
            None,
            Some(""),
        ];
        assert_eq!(most_common(v.into_iter()), "a");
        assert_eq!(most_common([Some("."), None].into_iter()), ".");
    }

    #[test]
    fn exact_int_compare() {
        assert!(float_ge_int(100_000.0, 100_000));
        assert!(!float_ge_int(99_999.5, 100_000));
        assert!(!float_ge_int(
            9_007_199_254_740_992.0,
            9_007_199_254_740_993
        ));
        assert!(!float_ge_int(f64::NAN, 0));
    }

    #[test]
    fn clean_cells() {
        assert_eq!(clean("  a\t b  "), "a b");
        assert_eq!(clean(" \t "), ".");
    }
}
