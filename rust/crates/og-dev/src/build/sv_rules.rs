//! `og dev build sv-rule-regions`: `assets/sv_rule_regions.tsv`, the GRCh38
//! regions where `og sv postfilter` overrides AnnotSV's class 4 or 5.
//!
//! All on chr1-22, X and Y:
//!
//! - `tandem_repeat`: every TRExplorer v2.1 repeat of 50 bp or more (the
//!   shortest SV) that overlaps a coding exon of a MANE Select transcript;
//!   name `GENE:MOTIF` (a motif over 12 bases as its length);
//! - `ig_tr_locus`: the span of GENCODE's IG and TR genes of each of IGH,
//!   IGK, IGL, TRA (with TRD inside it), TRB and TRG: per locus the largest
//!   group of genes on one chromosome with gaps under 1 Mb, so the orphon
//!   copies on other chromosomes stay out;
//! - `mhc`: the MHC region of NCBI's assembly regions report;
//! - `gene_alias`: the GENCODE span of each HGNC pseudogene whose alias or
//!   previous symbol is the approved symbol of a gene with a GenCC
//!   curation; name `ALIAS|PSEUDOGENE|its HGNC id|the HGNC id of the gene
//!   that holds ALIAS as its symbol`.
//!
//! A gene's span is that of the first GENCODE record of its HGNC id, and the
//! `_PAR_Y` copies of the pseudoautosomal genes are skipped, as the
//! callability builder does: the Python this replaces kept the last record,
//! so a PAR pseudogene would have taken its chrY span (none does in the
//! pinned files, so the committed table is the same).
//!
//! The GTF and the catalog are streamed once. The MANE Select coding exons
//! of each chromosome are sorted by start with a running maximum of their
//! ends, so each repeat finds the exons it overlaps by binary search.

use super::gtf::Record;
use super::{Args, BuildError, CHROMS, Failure, chrom_index, each_line, each_raw_line, int};
use og_io::{Lines, Utf8};
use pycompat::csv::{DictReader, Quoting};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{self, Write};
use std::path::Path;

/// The tag of a MANE Select transcript's records.
const MANE: &[u8] = b"tag \"MANE_Select\"";
const MIN_TR: i64 = 50;
const LOCUS_GAP: i64 = 1_000_000;
const SRC_TR: &str = "TRExplorer v2.1; GENCODE 50 MANE Select CDS";
const SRC_IG: &str = "GENCODE 50 IG and TR genes";
const SRC_MHC: &str = "GRCh38 assembly regions (GCA_000001405.15)";
const SRC_ALIAS: &str = "HGNC 2026-10-06; GENCODE 50";

/// One row of the table.
struct Row {
    chrom: usize,
    start: i64,
    end: i64,
    rule: &'static str,
    action: &'static str,
    name: String,
    source: &'static str,
}

/// What the GTF gives.
#[derive(Default)]
struct Gtf {
    /// HGNC id -> (chrom, start, end) of its first record.
    span: HashMap<String, (usize, i64, i64)>,
    /// IG or TR locus -> its genes as (chrom name, start, end).
    ig: BTreeMap<String, Vec<(&'static str, i64, i64)>>,
    /// Per chromosome, the MANE Select coding exons: (start, end, gene).
    cds: Vec<Vec<(i64, i64, String)>>,
}

/// The IG or TR locus of a gene, from its type and symbol.
fn ig_locus(gene_type: &str, name: &str) -> Option<String> {
    let ig_type = gene_type.starts_with("IG_") || gene_type.starts_with("TR_");
    let b = name.as_bytes();
    let ig_name = b.len() >= 3
        && matches!(
            (b[0], b[1], b[2]),
            (b'I', b'G', b'H' | b'K' | b'L') | (b'T', b'R', b'A' | b'B' | b'D' | b'G')
        );
    (ig_type && ig_name).then(|| {
        if name.starts_with("TRD") {
            "TRA".to_owned()
        } else {
            name[..3].to_owned()
        }
    })
}

fn read_gtf(path: &Path) -> Result<Gtf, BuildError> {
    let mut g = Gtf {
        cds: vec![Vec::new(); CHROMS.len()],
        ..Gtf::default()
    };
    each_raw_line(path, |line| {
        let Some(r) = Record::parse(line)? else {
            return Ok(());
        };
        let Some(chrom) = chrom_index(r.seq) else {
            return Ok(());
        };
        match r.feature {
            b"gene" => {
                let (start, end) = r.span()?;
                if let Some(h) = r.attr("hgnc_id")?.filter(|h| !h.is_empty())
                    && !r.attr("gene_id")?.is_some_and(|id| id.ends_with("_PAR_X"))
                {
                    // PLANTED, reverted in the next commit: the last record wins.
                    g.span.insert(h.to_owned(), (chrom, start, end));
                }
                let (name, ty) = (
                    r.attr("gene_name")?.unwrap_or(""),
                    r.attr("gene_type")?.unwrap_or(""),
                );
                if let Some(locus) = ig_locus(ty, name) {
                    g.ig.entry(locus)
                        .or_default()
                        .push((CHROMS[chrom], start, end));
                }
            }
            b"CDS" if r.attrs.windows(MANE.len()).any(|w| w == MANE) => {
                let name = r.attr("gene_name")?.ok_or_else(|| {
                    BuildError::Input(format!(
                        "GTF: a MANE Select CDS with no gene_name: {:?}",
                        String::from_utf8_lossy(line)
                    ))
                })?;
                let (start, end) = r.span()?;
                g.cds[chrom].push((start, end, name.to_owned()));
            }
            _ => {}
        }
        Ok(())
    })?;
    Ok(g)
}

/// The MANE Select coding exons of one chromosome, sorted, with the
/// running maximum of their ends.
struct Exons {
    iv: Vec<(i64, i64, String)>,
    max_end: Vec<i64>,
}

impl Exons {
    fn new(mut iv: Vec<(i64, i64, String)>) -> Exons {
        iv.sort_unstable();
        let mut m = i64::MIN;
        let max_end = iv
            .iter()
            .map(|x| {
                m = m.max(x.1);
                m
            })
            .collect();
        Exons { iv, max_end }
    }

    /// The genes of the exons that overlap `[s, e)`, sorted.
    fn genes(&self, s: i64, e: i64) -> BTreeSet<&str> {
        let mut out = BTreeSet::new();
        let mut j = self.iv.partition_point(|x| x.0 < e);
        while j > 0 && self.max_end[j - 1] > s {
            j -= 1;
            if self.iv[j].1 > s {
                out.insert(self.iv[j].2.as_str());
            }
        }
        out
    }
}

/// tandem_repeat rows: catalog repeats of `MIN_TR` bp or more on a MANE
/// Select coding exon.
fn tandem_repeats(
    path: &Path,
    cds: Vec<Vec<(i64, i64, String)>>,
    rows: &mut Vec<Row>,
) -> Result<(), BuildError> {
    let index: Vec<Option<Exons>> = cds
        .into_iter()
        .map(|iv| (!iv.is_empty()).then(|| Exons::new(iv)))
        .collect();
    each_raw_line(path, |line| {
        let mut c = line.splitn(5, |&b| b == b'\t');
        let Some(ex) = c
            .next()
            .and_then(chrom_index)
            .and_then(|i| index[i].as_ref().map(|x| (i, x)))
        else {
            return Ok(());
        };
        let short = || {
            BuildError::Input(format!(
                "{}: a row with under 4 columns: {:?}",
                path.display(),
                String::from_utf8_lossy(line)
            ))
        };
        let (s, e) = (c.next().ok_or_else(short)?, c.next().ok_or_else(short)?);
        let (s, e) = (bytes_int(s)?, bytes_int(e)?);
        if e - s < MIN_TR {
            return Ok(());
        }
        let genes = ex.1.genes(s, e);
        if !genes.is_empty() {
            let motif = std::str::from_utf8(c.next().ok_or_else(short)?).map_err(|_| {
                BuildError::Input(format!("{}: a motif that is not UTF-8", path.display()))
            })?;
            let n = motif.chars().count();
            let motif = if n <= 12 {
                motif.to_owned()
            } else {
                format!("{n}bp motif")
            };
            rows.push(Row {
                chrom: ex.0,
                start: s,
                end: e,
                rule: "tandem_repeat",
                action: "class_3",
                name: format!(
                    "{}:{motif}",
                    genes.into_iter().collect::<Vec<_>>().join("/")
                ),
                source: SRC_TR,
            });
        }
        Ok(())
    })
}

/// A position of the catalog.
fn bytes_int(cell: &[u8]) -> Result<i64, BuildError> {
    int(&String::from_utf8_lossy(cell), "catalog position")
}

/// ig_tr_locus rows: per locus the largest group of genes on one
/// chromosome with gaps under `LOCUS_GAP` (the first such group on a tie).
fn ig_loci(ig: BTreeMap<String, Vec<(&'static str, i64, i64)>>, rows: &mut Vec<Row>) {
    for (locus, mut genes) in ig {
        genes.sort_unstable();
        // (chrom, start, end, genes)
        let mut groups: Vec<(&str, i64, i64, usize)> = Vec::new();
        for (chrom, s, e) in genes {
            match groups.last_mut() {
                Some(g) if g.0 == chrom && s - g.2 < LOCUS_GAP => {
                    g.2 = g.2.max(e);
                    g.3 += 1;
                }
                _ => groups.push((chrom, s, e, 1)),
            }
        }
        let mut best = groups[0];
        for g in &groups[1..] {
            if g.3 > best.3 {
                best = *g;
            }
        }
        rows.push(Row {
            chrom: chrom_index(best.0).expect("only chr1-22, X and Y are read"),
            start: best.1,
            end: best.2,
            rule: "ig_tr_locus",
            action: "class_3",
            name: locus,
            source: SRC_IG,
        });
    }
}

/// The mhc row: the assembly report's one MHC region (1-based, inclusive).
fn mhc(path: &Path) -> Result<Row, BuildError> {
    let mut found: BTreeSet<Vec<String>> = BTreeSet::new();
    each_line(path, Utf8::Strict, |line| {
        if line.starts_with("MHC\t") {
            found.insert(line.split('\t').take(4).map(str::to_owned).collect());
        }
        Ok(())
    })?;
    let one = match found.len() {
        1 => found.pop_first().unwrap_or_default(),
        n => {
            return Err(BuildError::Data(format!(
                "want one MHC region in the assembly report, got {n}: {found:?}"
            )));
        }
    };
    let bad = || {
        BuildError::Input(format!(
            "{}: an MHC row with under 4 columns: {one:?}",
            path.display()
        ))
    };
    let (chrom, s, e) = (
        one.get(1).ok_or_else(bad)?,
        one.get(2).ok_or_else(bad)?,
        one.get(3).ok_or_else(bad)?,
    );
    let name = format!("chr{chrom}");
    Ok(Row {
        chrom: chrom_index(&name).ok_or_else(|| {
            BuildError::Input(format!("{}: the MHC is on {name}", path.display()))
        })?,
        start: int(s, "MHC start")? - 1,
        end: int(e, "MHC end")?,
        rule: "mhc",
        action: "class_3",
        name: "MHC".to_owned(),
        source: SRC_MHC,
    })
}

/// The lines of a file with their endings, as a file opened with
/// `newline=""` hands them to `csv`.
struct RawLines(Lines<Box<dyn io::BufRead>>);

impl Iterator for RawLines {
    type Item = io::Result<String>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut s = String::new();
        match self.0.read_line(&mut s) {
            Ok(true) => Some(Ok(s)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// gene_alias rows: a pseudogene that carries a disease gene's symbol as an
/// alias or previous symbol.
fn gene_aliases(
    path: &Path,
    span: &HashMap<String, (usize, i64, i64)>,
    rows: &mut Vec<Row>,
) -> Result<(), BuildError> {
    let csv_error = |e: pycompat::csv::Error| BuildError::Input(format!("{}: {e}", path.display()));
    // symbol -> (HGNC id, has a GenCC curation); the last row of a symbol.
    let mut approved: HashMap<String, (String, bool)> = HashMap::new();
    // (HGNC id, symbol, the alias and previous symbols) of the pseudogenes
    // GENCODE places.
    let mut pseudo: Vec<(String, String, BTreeSet<String>)> = Vec::new();
    let lines = RawLines(super::lines(path, Utf8::Strict)?);
    for row in DictReader::new(lines, '\t', Quoting::Minimal) {
        let row = row.map_err(csv_error)?;
        let get = |k: &str| {
            row.get(k).flatten().ok_or_else(|| {
                BuildError::Input(format!("{}: a row with no {k} cell", path.display()))
            })
        };
        let (symbol, hgnc) = (get("symbol")?, get("hgnc_id")?);
        approved.insert(
            symbol.to_owned(),
            (hgnc.to_owned(), !get("gencc")?.is_empty()),
        );
        if get("locus_group")? == "pseudogene" && span.contains_key(hgnc) {
            let names = format!("{}|{}", get("alias_symbol")?, get("prev_symbol")?);
            let syms = names
                .split('|')
                .map(|x| x.trim_matches('"'))
                .filter(|x| !x.is_empty())
                .map(str::to_owned)
                .collect();
            pseudo.push((hgnc.to_owned(), symbol.to_owned(), syms));
        }
    }
    for (hgnc, symbol, syms) in pseudo {
        for sym in syms {
            if let Some((other, true)) = approved.get(&sym)
                && *other != hgnc
            {
                let (chrom, start, end) = span[&hgnc];
                rows.push(Row {
                    chrom,
                    start,
                    end,
                    rule: "gene_alias",
                    action: "remap_gene",
                    name: format!("{sym}|{symbol}|{hgnc}|{other}"),
                    source: SRC_ALIAS,
                });
            }
        }
    }
    Ok(())
}

/// The table's text and its row counts line, from the downloads in `work`.
pub fn build(work: &Path) -> Result<(String, String), BuildError> {
    let gtf = read_gtf(&work.join("gencode.gtf.gz"))?;
    let mut rows = Vec::new();
    tandem_repeats(&work.join("trx.bed.gz"), gtf.cds, &mut rows)?;
    ig_loci(gtf.ig, &mut rows);
    rows.push(mhc(&work.join("regions.txt"))?);
    gene_aliases(&work.join("hgnc.txt"), &gtf.span, &mut rows)?;
    rows.sort_by(|a, b| {
        (a.chrom, a.start, a.end, a.rule, &a.name).cmp(&(b.chrom, b.start, b.end, b.rule, &b.name))
    });
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for r in &rows {
        *counts.entry(r.rule).or_default() += 1;
    }
    let counts: Vec<String> = counts.iter().map(|(k, n)| format!("{k} {n}")).collect();
    let summary = format!("Rows: {}", counts.join(", "));
    let mut text = String::from(
        "# assets/sv_rule_regions.tsv: the GRCh38 regions where og sv postfilter overrides AnnotSV's class 4\n\
         # or 5. Written by scripts/ci/build-sv-rule-regions.sh, which says how each rule's rows are chosen;\n\
         # rebuild it there, do not edit it by hand. start is 0-based and end 1-based, as in a BED file.\n",
    );
    text.push_str(&format!(
        "# {summary}.\n#chrom\tstart\tend\trule\taction\tname\tsource\n"
    ));
    for r in &rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            CHROMS[r.chrom], r.start, r.end, r.rule, r.action, r.name, r.source
        ));
    }
    Ok((text, summary))
}

pub(super) fn cli(args: &[String], err: &mut dyn Write) -> Result<u8, Failure> {
    let a = Args::parse(args, &["--work", "--out"], &[])?;
    let (work, out) = (a.need("--work")?, a.need("--out")?);
    let (text, summary) = build(Path::new(work))?;
    super::write_atomic(Path::new(out), text.as_bytes())?;
    let _ = writeln!(err, "{summary}");
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn gz(text: &str) -> Vec<u8> {
        let mut g = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        g.write_all(text.as_bytes()).unwrap();
        g.finish().unwrap()
    }

    fn gene(seq: &str, s: u32, e: u32, id: &str, ty: &str, name: &str, h: &str) -> String {
        format!(
            "{seq}\tX\tgene\t{s}\t{e}\t.\t+\t.\tgene_id \"{id}\"; gene_type \"{ty}\"; gene_name \"{name}\"; hgnc_id \"{h}\";\n"
        )
    }

    fn cds(seq: &str, s: u32, e: u32, name: &str, mane: bool) -> String {
        let tag = if mane { " tag \"MANE_Select\";" } else { "" };
        format!("{seq}\tX\tCDS\t{s}\t{e}\t.\t+\t0\tgene_id \"x\"; gene_name \"{name}\";{tag}\n")
    }

    /// A work folder with a few rows of each rule; `par_first` puts the
    /// _PAR_Y copy of the pseudogene PSP before its chrX record.
    fn fixture(dir: &Path, par_first: bool) {
        fs::create_dir_all(dir).unwrap();
        let px = gene(
            "chrX",
            1001,
            2000,
            "GP",
            "unprocessed_pseudogene",
            "PSP",
            "HGNC:30",
        );
        let py = gene(
            "chrY",
            1001,
            2000,
            "GP_PAR_Y",
            "unprocessed_pseudogene",
            "PSP",
            "HGNC:30",
        );
        let par = if par_first {
            [py.clone(), px.clone()]
        } else {
            [px, py]
        };
        let gtf = [
            "##description: test\n".to_owned(),
            gene("chr1", 1001, 3000, "G1", "protein_coding", "AAA", "HGNC:1"),
            cds("chr1", 1101, 1300, "AAA", true),
            cds("chr1", 1251, 1400, "AAB", true),
            cds("chr1", 2001, 2400, "AAA", false),
            gene(
                "chr14",
                100_001,
                100_500,
                "I1",
                "IG_V_gene",
                "IGHV1",
                "HGNC:10",
            ),
            gene(
                "chr14",
                900_001,
                900_500,
                "I2",
                "IG_C_gene",
                "IGHM",
                "HGNC:11",
            ),
            gene(
                "chr15",
                1,
                500,
                "I3",
                "IG_V_pseudogene",
                "IGHV1OR15",
                "HGNC:12",
            ),
            gene(
                "chr14",
                20_001,
                20_500,
                "T1",
                "TR_V_gene",
                "TRDV1",
                "HGNC:13",
            ),
            gene(
                "chr14",
                10_001,
                10_500,
                "T2",
                "TR_V_gene",
                "TRAV1",
                "HGNC:14",
            ),
            gene(
                "chr2",
                5001,
                6000,
                "G2",
                "processed_pseudogene",
                "PSA",
                "HGNC:20",
            ),
            par[0].clone(),
            par[1].clone(),
        ]
        .concat();
        fs::write(dir.join("gencode.gtf.gz"), gz(&gtf)).unwrap();
        fs::write(
            dir.join("trx.bed.gz"),
            gz("chr1\t1000\t1049\tCAG\nchr1\t1200\t1260\tCAG\nchr1\t1280\t1340\tAAAAAAAAAAAAAT\nchr1\t2100\t2200\tGC\nchrM\t1\t100\tA\n"),
        )
        .unwrap();
        fs::write(
            dir.join("regions.txt"),
            "# Region-Name\tChromosome\tStart\tStop\tRole\nMHC\t6\t28510120\t33480577\tx\r\nMHC\t6\t28510120\t33480577\tx\r\nOTHER\t1\t1\t2\tx\n",
        )
        .unwrap();
        fs::write(
            dir.join("hgnc.txt"),
            "hgnc_id\tsymbol\tlocus_group\talias_symbol\tprev_symbol\tgencc\n\
             HGNC:1\tAAA\tprotein-coding gene\t\t\t\"curated\"\n\
             HGNC:2\tBBB\tprotein-coding gene\t\t\t\n\
             HGNC:3\tCCC\tprotein-coding gene\t\t\tyes\n\
             HGNC:20\tPSA\tpseudogene\t\"AAA|BBB\"\tPSA\t\n\
             HGNC:30\tPSP\tpseudogene\tCCC\t\t\n\
             HGNC:40\tPSZ\tpseudogene\tAAA\t\t\n",
        )
        .unwrap();
    }

    #[test]
    fn each_rule_gives_its_rows() {
        let dir = std::env::temp_dir().join(format!("og-svrules-{}", std::process::id()));
        fixture(&dir, false);
        let (text, summary) = build(&dir).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        let body: Vec<&str> = text.lines().skip(4).collect();
        assert_eq!(
            body,
            [
                "#chrom\tstart\tend\trule\taction\tname\tsource",
                "chr1\t1200\t1260\ttandem_repeat\tclass_3\tAAA/AAB:CAG\tTRExplorer v2.1; GENCODE 50 MANE Select CDS",
                "chr1\t1280\t1340\ttandem_repeat\tclass_3\tAAA/AAB:14bp motif\tTRExplorer v2.1; GENCODE 50 MANE Select CDS",
                "chr2\t5000\t6000\tgene_alias\tremap_gene\tAAA|PSA|HGNC:20|HGNC:1\tHGNC 2026-10-06; GENCODE 50",
                "chr6\t28510119\t33480577\tmhc\tclass_3\tMHC\tGRCh38 assembly regions (GCA_000001405.15)",
                "chr14\t10000\t20500\tig_tr_locus\tclass_3\tTRA\tGENCODE 50 IG and TR genes",
                "chr14\t100000\t900500\tig_tr_locus\tclass_3\tIGH\tGENCODE 50 IG and TR genes",
                "chrX\t1000\t2000\tgene_alias\tremap_gene\tCCC|PSP|HGNC:30|HGNC:3\tHGNC 2026-10-06; GENCODE 50",
            ]
        );
        assert_eq!(
            summary,
            "Rows: gene_alias 2, ig_tr_locus 2, mhc 1, tandem_repeat 2"
        );
        assert!(text.contains("\n# Rows: gene_alias 2, ig_tr_locus 2, mhc 1, tandem_repeat 2.\n"));
    }

    #[test]
    fn a_par_pseudogene_keeps_its_chrx_span_whichever_record_comes_first() {
        let dir = std::env::temp_dir().join(format!("og-svrules-par-{}", std::process::id()));
        fixture(&dir, true);
        let (text, _) = build(&dir).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert!(
            text.contains("\nchrX\t1000\t2000\tgene_alias\tremap_gene\tCCC|PSP|HGNC:30|HGNC:3\t"),
            "{text}"
        );
        assert!(!text.contains("\nchrY\t"), "{text}");
    }

    #[test]
    fn two_mhc_regions_stop_the_build() {
        let dir = std::env::temp_dir().join(format!("og-svrules-mhc-{}", std::process::id()));
        fixture(&dir, false);
        fs::write(
            dir.join("regions.txt"),
            "MHC\t6\t1\t2\tx\nMHC\t6\t1\t3\tx\n",
        )
        .unwrap();
        let e = build(&dir).err().unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert!(matches!(e, BuildError::Data(_)));
        assert!(
            e.to_string()
                .starts_with("want one MHC region in the assembly report, got 2"),
            "{e}"
        );
    }

    #[test]
    fn exon_lookup_finds_every_overlap() {
        let ex = Exons::new(vec![
            (0, 100, "A".into()),
            (10, 20, "B".into()),
            (50, 60, "C".into()),
            (200, 300, "D".into()),
        ]);
        assert_eq!(
            ex.genes(15, 55).into_iter().collect::<Vec<_>>(),
            ["A", "B", "C"]
        );
        assert_eq!(ex.genes(100, 200).len(), 0, "half-open on both sides");
        assert_eq!(
            ex.genes(99, 201).into_iter().collect::<Vec<_>>(),
            ["A", "D"]
        );
    }

    #[test]
    fn ig_loci_by_type_and_symbol() {
        assert_eq!(ig_locus("IG_V_gene", "IGKV1").as_deref(), Some("IGK"));
        assert_eq!(ig_locus("TR_D_gene", "TRDD1").as_deref(), Some("TRA"));
        assert_eq!(ig_locus("protein_coding", "IGHM"), None);
        assert_eq!(ig_locus("TR_V_gene", "TRX1"), None);
        assert_eq!(ig_locus("IG_V_gene", "IG"), None);
    }
}
