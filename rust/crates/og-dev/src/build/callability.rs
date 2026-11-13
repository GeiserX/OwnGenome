//! `og dev build callability`: `assets/callability_genes.tsv`, the genes
//! where short reads are often not placed with confidence, so a "no
//! finding" there needs a reading of its own (step 16b measures each one).
//!
//! A gene is a GENCODE 50 gene with an HGNC id, on chr1-22, X or Y unless it
//! is ALT-only, and protein-coding unless it is a CMRG gene. It is listed
//! with every reason that applies:
//!
//! - `CMRG`: one of GIAB's 273 challenging medically relevant genes,
//!   matched by symbol, or by the protein-coding gene that covers most of
//!   its CMRG span when GENCODE has renamed it;
//! - `segdup`: half or more of the gene body in a segmental duplication;
//! - `false_duplication`: the gene body overlaps a false duplication of
//!   GRCh38 (either copy);
//! - `collapsed`: the gene body overlaps a collapsed duplication;
//! - `ALT-only`: the gene is on an ALT scaffold and nowhere on chr1-22, X,
//!   Y or M;
//!
//! and the step that calls it with a method of its own: `parascopy` for the
//! genes of `assets/paralog_loci.tsv`, `pypgx` for the bam genes of
//! `assets/pgx_genes.tsv` callable without ALT contigs, `cn_loci:<locus>`
//! for a gene a gene row of `bin/cn_loci.tsv` names; `.` for none.
//!
//! `--work` holds the seven downloads `scripts/ci/build-callability.sh`
//! names. The GTF is read once, line by line; each interval file is merged
//! once per chromosome, so a gene's coverage is a binary search.

use super::gtf::Record;
use super::{Args, BuildError, CHROMS, Failure, chrom_index, each_line, each_raw_line, int};
use og_intervals::sets::Set;
use og_io::Utf8;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

/// The reasons, in the order a row lists them.
const REASONS: [&str; 5] = [
    "CMRG",
    "segdup",
    "false_duplication",
    "collapsed",
    "ALT-only",
];
const CMRG: u8 = 1;
const SEGDUP: u8 = 1 << 1;
const FALSE_DUP: u8 = 1 << 2;
const COLLAPSED: u8 = 1 << 3;
const ALT_ONLY: u8 = 1 << 4;

/// A gene of GENCODE on chr1-22, X or Y: the first record of its HGNC id.
struct Gene {
    symbol: String,
    hgnc: String,
    chrom: usize,
    start: i64,
    end: i64,
    coding: bool,
}

/// What the GTF gives.
#[derive(Default)]
struct Genes {
    /// In the GTF's order.
    primary: Vec<Gene>,
    by_hgnc: HashMap<String, usize>,
    /// Every gene symbol on chr1-22, X, Y or M.
    on_primary: HashSet<String>,
    /// HGNC id -> (symbol, UCSC ALT scaffold name) of the first
    /// protein-coding ALT copy.
    alt: HashMap<String, (String, String)>,
}

/// The ALT scaffolds of NCBI's assembly report: GenBank accession -> UCSC
/// name (GENCODE names them by accession).
fn read_alt(path: &Path) -> Result<HashMap<String, String>, BuildError> {
    let mut alt = HashMap::new();
    each_line(path, Utf8::Strict, |line| {
        let c: Vec<&str> = line.split('\t').collect();
        if !line.starts_with('#') && c.len() >= 10 && c[1] == "alt-scaffold" {
            alt.insert(c[4].to_owned(), c[9].to_owned());
        }
        Ok(())
    })?;
    Ok(alt)
}

fn read_gtf(path: &Path, alt: &HashMap<String, String>) -> Result<Genes, BuildError> {
    let mut g = Genes::default();
    each_raw_line(path, |line| {
        let Some(r) = Record::parse(line)? else {
            return Ok(());
        };
        if r.feature != b"gene" {
            return Ok(());
        }
        let chrom = chrom_index(r.seq);
        if (chrom.is_some() || r.seq == b"chrM")
            && let Some(n) = r.attr("gene_name")?
        {
            g.on_primary.insert(n.to_owned());
        }
        let Some(h) = r.attr("hgnc_id")?.filter(|h| !h.is_empty()) else {
            return Ok(());
        };
        let need = |k: &str| {
            r.attr(k)?.ok_or_else(|| {
                BuildError::Input(format!(
                    "GTF: a gene with an HGNC id and no {k}: {:?}",
                    String::from_utf8_lossy(line)
                ))
            })
        };
        if need("gene_id")?.ends_with("_PAR_Y") {
            return Ok(());
        }
        let coding = r.attr("gene_type")? == Some("protein_coding");
        match chrom {
            Some(chrom) if !g.by_hgnc.contains_key(h) => {
                let (start, end) = r.span()?;
                g.by_hgnc.insert(h.to_owned(), g.primary.len());
                g.primary.push(Gene {
                    symbol: need("gene_name")?.to_owned(),
                    hgnc: h.to_owned(),
                    chrom,
                    start,
                    end,
                    coding,
                });
            }
            _ => {
                if let Some(name) = std::str::from_utf8(r.seq).ok().and_then(|s| alt.get(s))
                    && coding
                    && !g.alt.contains_key(h)
                {
                    g.alt
                        .insert(h.to_owned(), (need("gene_name")?.to_owned(), name.clone()));
                }
            }
        }
        Ok(())
    })?;
    Ok(g)
}

/// A BED file as a [`Set`] per sequence.
fn read_sets(paths: &[&Path]) -> Result<HashMap<String, Set>, BuildError> {
    let mut rows: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for path in paths {
        each_line(path, Utf8::Strict, |line| {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() >= 3
                && !(line.starts_with('#')
                    || line.starts_with("track")
                    || line.starts_with("browser"))
            {
                let s = u64::try_from(int(c[1], "BED start")?);
                let e = u64::try_from(int(c[2], "BED end")?);
                let (Ok(s), Ok(e)) = (s, e) else {
                    return Err(BuildError::Input(format!(
                        "{}: a negative position: {line:?}",
                        path.display()
                    )));
                };
                rows.entry(c[0].to_owned()).or_default().push((s, e));
            }
            Ok(())
        })?;
    }
    Ok(rows.into_iter().map(|(k, v)| (k, Set::new(v))).collect())
}

/// Bases of `[s, e)` on `chrom` inside `sets`.
fn covered(sets: &HashMap<String, Set>, chrom: usize, s: i64, e: i64) -> u64 {
    let (s, e) = (u64::try_from(s).unwrap_or(0), u64::try_from(e).unwrap_or(0));
    sets.get(CHROMS[chrom]).map_or(0, |set| set.covered(s, e))
}

/// The rows of a tab-separated table after its comment lines, as
/// `csv.DictReader` with `QUOTE_NONE` reads them: the first line is the
/// header; each row is handed to `f` with a cell getter.
fn each_row(
    path: &Path,
    mut f: impl FnMut(&dyn Fn(&str) -> Result<String, BuildError>) -> Result<(), BuildError>,
) -> Result<(), BuildError> {
    let mut header: Option<Vec<String>> = None;
    each_line(path, Utf8::Strict, |line| {
        if pycompat::text::strip(line, None).is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let Some(h) = &header else {
            header = Some(cells.iter().map(|c| (*c).to_owned()).collect());
            return Ok(());
        };
        let get = |k: &str| {
            h.iter()
                .rposition(|n| n == k)
                .and_then(|i| cells.get(i))
                .map(|c| (*c).to_owned())
                .ok_or_else(|| {
                    BuildError::Input(format!("{}: no {k} cell in {line:?}", path.display()))
                })
        };
        f(&get)
    })
}

/// Who calls a gene with a method of its own.
struct Owners {
    parascopy: HashSet<String>,
    pypgx: HashSet<String>,
    /// A name a `bin/cn_loci.tsv` gene row reports -> its locus (the first
    /// row that names it).
    cn: HashMap<String, String>,
}

impl Owners {
    fn read(repo: &Path) -> Result<Owners, BuildError> {
        let mut parascopy = HashSet::new();
        each_row(&repo.join("assets/paralog_loci.tsv"), |get| {
            parascopy.extend(get("genes")?.split(',').map(str::to_owned));
            Ok(())
        })?;
        // A gene step 32 skips on the default no-ALT reference
        // (callable_without_alt no) has no caller there, and step 32
        // refuses a BAM with ALT contigs, so it is not pypgx's.
        let mut pypgx = HashSet::new();
        each_row(&repo.join("assets/pgx_genes.tsv"), |get| {
            if get("mode")? == "bam" && get("callable_without_alt")? != "no" {
                pypgx.insert(get("gene")?);
            }
            Ok(())
        })?;
        // A gene row reports the genes it names (its locus, as RHD, or LPA
        // for LPA_KIV2, or a target, as HBA1 and HBA2), not every gene whose
        // span holds one of its targets.
        let mut cn: HashMap<String, String> = HashMap::new();
        each_row(&repo.join("bin/cn_loci.tsv"), |get| {
            if get("kind")? != "gene" {
                return Ok(());
            }
            let locus = get("locus")?;
            let targets = get("targets")?;
            let mut names = vec![locus.as_str(), first(&locus, '_')];
            for part in targets.split(';') {
                let n = first(part, ':');
                names.extend([n, first(n, '_')]);
            }
            for n in names {
                cn.entry(n.to_owned()).or_insert_with(|| locus.clone());
            }
            Ok(())
        })?;
        Ok(Owners {
            parascopy,
            pypgx,
            cn,
        })
    }

    fn of(&self, symbol: &str) -> String {
        if self.parascopy.contains(symbol) {
            "parascopy".to_owned()
        } else if self.pypgx.contains(symbol) {
            "pypgx".to_owned()
        } else if let Some(locus) = self.cn.get(symbol) {
            format!("cn_loci:{locus}")
        } else {
            ".".to_owned()
        }
    }
}

/// The text before the first `sep`.
fn first(s: &str, sep: char) -> &str {
    s.split(sep).next().unwrap_or(s)
}

/// The reasons of each gene, by HGNC id.
fn reasons(work: &Path, genes: &Genes) -> Result<HashMap<String, u8>, BuildError> {
    let mut why: HashMap<String, u8> = HashMap::new();
    let mut by_symbol: HashMap<&str, usize> = HashMap::new();
    for (i, g) in genes.primary.iter().enumerate() {
        by_symbol.insert(&g.symbol, i);
    }
    let cmrg = work.join("cmrg.bed");
    each_line(&cmrg, Utf8::Strict, |line| {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 4 {
            return Err(BuildError::Input(format!(
                "cmrg.bed: a row with under 4 columns: {line:?}"
            )));
        }
        let (chrom, sym) = (c[0], c[3]);
        let (s, e) = (int(c[1], "cmrg.bed start")?, int(c[2], "cmrg.bed end")?);
        let i = match by_symbol.get(sym) {
            Some(&i) => i,
            None => {
                // Renamed since CMRG v1.00: the protein-coding gene that
                // covers most of its span (the larger HGNC id on a tie).
                let best = genes
                    .primary
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| g.coding && CHROMS[g.chrom] == chrom)
                    .map(|(i, g)| (e.min(g.end) - s.max(g.start), &g.hgnc, i))
                    .max_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
                match best {
                    Some((overlap, _, i)) if overlap > 0 => i,
                    _ => {
                        return Err(BuildError::Data(format!(
                            "CMRG gene {sym} ({chrom}:{s}-{e}) matches no GENCODE 50 gene"
                        )));
                    }
                }
            }
        };
        *why.entry(genes.primary[i].hgnc.clone()).or_default() |= CMRG;
        Ok(())
    })?;
    let segdups = read_sets(&[&work.join("segdups.bed.gz")])?;
    let collapsed = read_sets(&[&work.join("collapsed.bed.gz")])?;
    let falsedup = read_sets(&[
        &work.join("falsedup_correct.bed.gz"),
        &work.join("falsedup_incorrect.bed.gz"),
    ])?;
    for g in genes.primary.iter().filter(|g| g.coding) {
        let mut bits = 0;
        // Half or more: covered >= 0.5 * length, in whole numbers.
        if 2 * covered(&segdups, g.chrom, g.start, g.end)
            >= u64::try_from(g.end - g.start).unwrap_or(0)
        {
            bits |= SEGDUP;
        }
        if covered(&falsedup, g.chrom, g.start, g.end) > 0 {
            bits |= FALSE_DUP;
        }
        if covered(&collapsed, g.chrom, g.start, g.end) > 0 {
            bits |= COLLAPSED;
        }
        if bits != 0 {
            *why.entry(g.hgnc.clone()).or_default() |= bits;
        }
    }
    for (h, (sym, _)) in &genes.alt {
        if !genes.by_hgnc.contains_key(h) && !genes.on_primary.contains(sym) {
            *why.entry(h.clone()).or_default() |= ALT_ONLY;
        }
    }
    Ok(why)
}

/// A row of the catalog: its sort key (chrom order, start, symbol, HGNC
/// id; ALT-only genes last), its cells and its reasons.
type CatalogRow<'a> = ((usize, i64, &'a str, &'a str), [String; 7], u8);

/// The catalog's text and the summary line, from the downloads in `work`
/// and the owner tables of `repo`.
pub fn build(work: &Path, repo: &Path) -> Result<(String, String), BuildError> {
    let alt = read_alt(&work.join("report.txt"))?;
    let genes = read_gtf(&work.join("gencode.gtf.gz"), &alt)?;
    let why = reasons(work, &genes)?;
    let owners = Owners::read(repo)?;
    let mut rows: Vec<CatalogRow> = Vec::with_capacity(why.len());
    for (h, &bits) in &why {
        let (key, chrom, start, end, sym) = match genes.by_hgnc.get(h) {
            Some(&i) => {
                let g = &genes.primary[i];
                let key = (g.chrom, g.start, g.symbol.as_str(), h.as_str());
                (
                    key,
                    CHROMS[g.chrom].to_owned(),
                    g.start.to_string(),
                    g.end.to_string(),
                    &g.symbol,
                )
            }
            None => {
                let (sym, name) = &genes.alt[h];
                (
                    (CHROMS.len(), 0, sym.as_str(), h.as_str()),
                    name.clone(),
                    ".".to_owned(),
                    ".".to_owned(),
                    sym,
                )
            }
        };
        let listed: Vec<&str> = REASONS
            .iter()
            .enumerate()
            .filter(|(i, _)| bits & (1 << i) != 0)
            .map(|(_, r)| *r)
            .collect();
        let cells = [
            sym.clone(),
            h.clone(),
            chrom,
            start,
            end,
            listed.join(","),
            owners.of(sym),
        ];
        rows.push((key, cells, bits));
    }
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let mut seen: HashMap<&str, usize> = HashMap::new();
    for r in &rows {
        *seen.entry(r.0.2).or_default() += 1;
    }
    let mut dup: Vec<&str> = seen
        .into_iter()
        .filter(|&(_, n)| n > 1)
        .map(|(s, _)| s)
        .collect();
    if !dup.is_empty() {
        dup.sort_unstable();
        return Err(BuildError::Data(format!(
            "gene symbols listed twice: {}",
            dup.join(", ")
        )));
    }
    let counts: Vec<String> = REASONS
        .iter()
        .enumerate()
        .map(|(i, r)| {
            format!(
                "{r} {}",
                rows.iter().filter(|x| x.2 & (1 << i) != 0).count()
            )
        })
        .collect();
    let owned = rows.iter().filter(|r| r.1[6] != ".").count();
    let summary = format!(
        "Genes: {} ({}); with an owner: {owned}",
        rows.len(),
        counts.join(", ")
    );
    let mut text = String::from(
        "# assets/callability_genes.tsv: the genes where short reads are often not placed with confidence, so a\n\
         # 'no finding' there needs a reading of its own. Step 16b measures each gene on chr1-22, X and Y; the\n\
         # report lists the ones it could not read. Written by scripts/ci/build-callability.sh, which says how the\n\
         # rows and the owners are chosen; rebuild it there, do not edit it by hand. start is 0-based and end\n\
         # 1-based, as in a BED file; an ALT-only gene has its ALT scaffold as chrom and '.' as start and end.\n",
    );
    let _ = writeln!(text, "# {summary}.");
    text.push_str("gene\thgnc_id\tchrom\tstart\tend\treasons\towner\n");
    for r in &rows {
        text.push_str(&r.1.join("\t"));
        text.push('\n');
    }
    Ok((text, summary))
}

pub(super) fn cli(args: &[String], err: &mut dyn Write) -> Result<u8, Failure> {
    let a = Args::parse(args, &["--work", "--repo", "--out"], &[])?;
    let (work, repo, out) = (a.need("--work")?, a.need("--repo")?, a.need("--out")?);
    let (text, summary) = build(Path::new(work), Path::new(repo))?;
    super::write_atomic(Path::new(out), text.as_bytes())?;
    let _ = writeln!(err, "{summary}");
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// A work folder and a repo folder with tiny inputs: A is a CMRG gene by
    /// symbol, OLD a CMRG gene GENCODE renamed to B, C sits in a segmental
    /// duplication, D over a false and a collapsed duplication, E is only on
    /// an ALT scaffold, F is on an ALT scaffold and chr1, P is a CMRG gene whose
    /// _PAR_Y copy comes first and must not move it to chrY, N is not protein-coding and has no reason.
    fn fixture(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
        let (work, repo) = (dir.join("work"), dir.join("repo"));
        fs::create_dir_all(&work).unwrap();
        fs::create_dir_all(repo.join("assets")).unwrap();
        fs::create_dir_all(repo.join("bin")).unwrap();
        let gene = |seq: &str, s: u32, e: u32, id: &str, ty: &str, name: &str, h: &str| {
            let hg = if h.is_empty() {
                String::new()
            } else {
                format!(" hgnc_id \"{h}\";")
            };
            format!(
                "{seq}\tX\tgene\t{s}\t{e}\t.\t+\t.\tgene_id \"{id}\"; gene_type \"{ty}\"; gene_name \"{name}\";{hg}\n"
            )
        };
        let gtf = [
            "##format: gtf\n".to_owned(),
            gene("chr1", 101, 200, "G1", "protein_coding", "A", "HGNC:1"),
            "chr1\tX\texon\t101\t150\t.\t+\t.\tgene_id \"G1\";\n".to_owned(),
            gene("chr1", 1001, 2000, "G2", "protein_coding", "B", "HGNC:2"),
            gene("chr2", 101, 1100, "G3", "protein_coding", "C", "HGNC:3"),
            gene("chr2", 5001, 6000, "G4", "protein_coding", "D", "HGNC:4"),
            gene("chr2", 9001, 9100, "G5", "lncRNA", "N", "HGNC:5"),
            gene(
                "chrY",
                100_001,
                100_500,
                "G6_PAR_Y",
                "protein_coding",
                "P",
                "HGNC:6",
            ),
            gene(
                "chrX",
                100_001,
                100_500,
                "G6",
                "protein_coding",
                "P",
                "HGNC:6",
            ),
            gene("chr1", 3001, 3100, "G7", "protein_coding", "F", "HGNC:7"),
            gene("KI1.1", 1, 100, "G8", "protein_coding", "E", "HGNC:8"),
            gene("KI1.1", 1, 100, "G9", "protein_coding", "F", "HGNC:9"),
            gene("chrM", 1, 50, "G10", "protein_coding", "MT", ""),
        ]
        .concat();
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(gtf.as_bytes()).unwrap();
        fs::write(work.join("gencode.gtf.gz"), gz.finish().unwrap()).unwrap();
        fs::write(
            work.join("report.txt"),
            "# report\nKI1\talt-scaffold\t1\tc\tKI1.1\t=\tx\ty\tz\tchr1_KI1v1_alt\r\n",
        )
        .unwrap();
        fs::write(
            work.join("cmrg.bed"),
            "chr1\t100\t200\tA\nchr1\t900\t1500\tOLD\nchrX\t100000\t100500\tP\n",
        )
        .unwrap();
        fs::write(
            work.join("segdups.bed.gz"),
            "chr2\t100\t400\nchr2\t350\t700\n",
        )
        .unwrap();
        fs::write(work.join("falsedup_correct.bed.gz"), "chr2\t5999\t6001\n").unwrap();
        fs::write(work.join("falsedup_incorrect.bed.gz"), "").unwrap();
        fs::write(work.join("collapsed.bed.gz"), "#c\nchr2\t5000\t5001\n").unwrap();
        fs::write(
            repo.join("assets/paralog_loci.tsv"),
            "# x\nlocus\tgenes\nl\tC,Z\n",
        )
        .unwrap();
        fs::write(
            repo.join("assets/pgx_genes.tsv"),
            "gene\tmode\tcallable_without_alt\nD\tbam\tyes\nA\tbam\tno\n",
        )
        .unwrap();
        fs::write(
            repo.join("bin/cn_loci.tsv"),
            "locus\tkind\ttargets\nA_x\tgene\tq_t1:1-2\nY\tinterval\tA:1-2\n",
        )
        .unwrap();
        (work, repo)
    }

    #[test]
    fn the_tiny_catalog_has_every_reason_and_owner() {
        let dir = std::env::temp_dir().join(format!("og-callability-{}", std::process::id()));
        let (work, repo) = fixture(&dir);
        let (text, summary) = build(&work, &repo).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            body,
            [
                "gene\thgnc_id\tchrom\tstart\tend\treasons\towner",
                "A\tHGNC:1\tchr1\t100\t200\tCMRG\tcn_loci:A_x",
                "B\tHGNC:2\tchr1\t1000\t2000\tCMRG\t.",
                "C\tHGNC:3\tchr2\t100\t1100\tsegdup\tparascopy",
                "D\tHGNC:4\tchr2\t5000\t6000\tfalse_duplication,collapsed\tpypgx",
                "P\tHGNC:6\tchrX\t100000\t100500\tCMRG\t.",
                "E\tHGNC:8\tchr1_KI1v1_alt\t.\t.\tALT-only\t.",
            ]
        );
        assert_eq!(
            summary,
            "Genes: 6 (CMRG 3, segdup 1, false_duplication 1, collapsed 1, ALT-only 1); with an owner: 3"
        );
        assert!(text.contains(&format!("\n# {summary}.\n")));
    }

    #[test]
    fn a_cmrg_gene_with_no_gencode_gene_stops_the_build() {
        let dir = std::env::temp_dir().join(format!("og-callability-cmrg-{}", std::process::id()));
        let (work, repo) = fixture(&dir);
        fs::write(work.join("cmrg.bed"), "chr3\t1\t2\tGONE\n").unwrap();
        let e = build(&work, &repo).err().unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            e.to_string(),
            "CMRG gene GONE (chr3:1-2) matches no GENCODE 50 gene"
        );
        assert!(matches!(e, BuildError::Data(_)));
    }

    #[test]
    fn a_symbol_listed_twice_stops_the_build() {
        let dir = std::env::temp_dir().join(format!("og-callability-dup-{}", std::process::id()));
        let (work, repo) = fixture(&dir);
        // Two genes named A, both CMRG (the second by its span).
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(b"chr1\tX\tgene\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_type \"protein_coding\"; gene_name \"A\"; hgnc_id \"HGNC:1\";\nchr1\tX\tgene\t901\t1500\t.\t+\t.\tgene_id \"G2\"; gene_type \"protein_coding\"; gene_name \"A\"; hgnc_id \"HGNC:2\";\n").unwrap();
        fs::write(work.join("gencode.gtf.gz"), gz.finish().unwrap()).unwrap();
        fs::write(
            work.join("cmrg.bed"),
            "chr1\t100\t200\tX1\nchr1\t900\t1500\tX2\n",
        )
        .unwrap();
        let e = build(&work, &repo).err().unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert_eq!(e.to_string(), "gene symbols listed twice: A");
    }

    #[test]
    fn coverage_is_union_coverage() {
        let sets = HashMap::from([(
            "chr2".to_owned(),
            Set::new([(100, 400), (350, 700), (900, 950)]),
        )]);
        assert_eq!(covered(&sets, 1, 100, 1100), 650);
        assert_eq!(covered(&sets, 0, 100, 1100), 0);
    }
}
