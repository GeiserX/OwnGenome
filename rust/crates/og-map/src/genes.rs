//! The map's gene files (`og map genes`): every gene's name and span, and
//! the exon blocks, per chromosome, built once per GENCODE release.
//!
//! - Inputs: setup's `reference/gencode.v<N>.basic.genes.gtf` (the gene
//!   lines of GENCODE's basic annotation) and `reference/gencode.v<N>.basic.exons.bed`
//!   (every exon of a basic transcript, once each, with no name), under
//!   GENOME_DIR. When several releases are installed, the highest `N` is read.
//! - GTF positions are 1-based and inclusive, so `s = start - 1` and
//!   `e = end`. The BED is already 0-based with an exclusive end.
//! - The exons become their union per chromosome: overlapping or touching
//!   exons merge into one block. A gene draws the blocks inside its span;
//!   a block is "an exon of some basic transcript here", not one transcript.
//! - Only the map's 24 chromosomes are kept (chrM and other contigs are
//!   skipped), so a name found in the names file always lands on a row.
//! - The files are cached under `<cache-dir>/genes/<stamp>/`, where the
//!   stamp covers both inputs' path, size and modification time, the format
//!   version and the og tree hash. They are written to a hidden folder and
//!   renamed into place in one step, so a reader never sees half a set.
//!   Writing a new set removes the older ones.
//!
//! `genes.<chrom>.json.gz`, one per chromosome, genes and blocks sorted by
//! start, starts as gaps from the previous start and ends as lengths:
//!
//! ```text
//! {"v":1,"gencode":"v50","chrom":"chr7",
//!  "g":{"ds":[],"l":[],"strand":[],"name":[],"type":[]},
//!  "x":{"ds":[],"l":[]}}
//! ```
//!
//! `genes.names.json.gz`, every gene sorted by name (byte order), then by
//! chromosome and start, with absolute spans:
//!
//! ```text
//! {"v":1,"gencode":"v50","name":[],"chrom":[],"s":[],"e":[]}
//! ```

use crate::doc::{gaps, obj};
use crate::reference::{CHROMS, chrom_index};
use crate::{BuildError, VERSION};
use flate2::{Compression, GzBuilder};
use pycompat::json::{self, Dump, Value};
use std::fs;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

/// What a reader is told when the GENCODE files are not installed.
pub const MISSING: &str = "Gene names need step 08's GENCODE file";

/// The file of the names, beside the per-chromosome files.
pub const NAMES_FILE: &str = "genes.names.json.gz";

/// Which file a request wants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Want {
    /// Index into [`CHROMS`].
    Chrom(usize),
    Names,
}

impl Want {
    /// `chrom` as one of the map's chromosomes, or a usage error.
    pub fn chrom(name: &str) -> Result<Want, BuildError> {
        chrom_index(name).map(Want::Chrom).ok_or_else(|| {
            BuildError::Usage(format!(
                "not a chromosome of the map: {name:?} (chr1-chr22, chrX or chrY)"
            ))
        })
    }

    pub fn file_name(self) -> String {
        match self {
            Want::Chrom(c) => format!("genes.{}.json.gz", CHROMS[c]),
            Want::Names => NAMES_FILE.to_owned(),
        }
    }
}

/// One gene line of the GTF.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    /// Index into [`CHROMS`].
    pub chrom: usize,
    /// 0-based start.
    pub s: u64,
    /// Exclusive end.
    pub e: u64,
    /// `+`, `-` or `.`.
    pub strand: u8,
    /// `gene_name`, else `gene_id`.
    pub name: String,
    /// `gene_type`, empty when the line has none.
    pub kind: String,
}

/// One release's genes and exon blocks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genes {
    /// `v50` for `gencode.v50.basic.genes.gtf`.
    pub release: String,
    /// Sorted by chromosome, start, end and name.
    pub genes: Vec<Gene>,
    /// Per chromosome in [`CHROMS`] order: the merged blocks, sorted.
    pub exons: Vec<Vec<(u64, u64)>>,
}

/// The two input files of one release.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Inputs {
    pub gtf: PathBuf,
    pub bed: PathBuf,
    pub release: String,
}

/// The release number of a genes GTF file name, if it is one.
fn release_of(file_name: &str) -> Option<u32> {
    let n = file_name
        .strip_prefix("gencode.v")?
        .strip_suffix(".basic.genes.gtf")?;
    if n.is_empty() || !n.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    n.parse().ok()
}

/// Finds the highest installed release under `<root>/reference`.
pub fn locate(root: &Path) -> Result<Inputs, BuildError> {
    let dir = root.join("reference");
    let newest = fs::read_dir(&dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter(|e| e.path().is_file())
        .filter_map(|e| release_of(&e.file_name().to_string_lossy()))
        .max();
    let Some(n) = newest else {
        return Err(BuildError::NoMap(format!(
            "{MISSING}: no reference/gencode.v<N>.basic.genes.gtf under the root (scripts/setup.sh installs it)"
        )));
    };
    let bed = dir.join(format!("gencode.v{n}.basic.exons.bed"));
    if !bed.is_file() {
        return Err(BuildError::NoMap(format!(
            "{MISSING}: reference/gencode.v{n}.basic.exons.bed is missing (scripts/setup.sh installs it)"
        )));
    }
    Ok(Inputs {
        gtf: dir.join(format!("gencode.v{n}.basic.genes.gtf")),
        bed,
        release: format!("v{n}"),
    })
}

/// Reads `input` line by line, handing each line without its end to `row`
/// with its 1-based number. Errors read `name:line: reason`.
fn each_line(
    mut input: impl BufRead,
    name: &str,
    mut row: impl FnMut(&str) -> Result<(), String>,
) -> Result<(), String> {
    let mut buf = Vec::with_capacity(512);
    let mut line_no = 0usize;
    loop {
        buf.clear();
        let n = input
            .read_until(b'\n', &mut buf)
            .map_err(|e| format!("{name}:{}: {e}", line_no + 1))?;
        if n == 0 {
            return Ok(());
        }
        line_no += 1;
        let mut line = &buf[..];
        while let [rest @ .., b'\n' | b'\r'] = line {
            line = rest;
        }
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let text = std::str::from_utf8(line).map_err(|_| format!("{name}:{line_no}: not UTF-8"))?;
        row(text).map_err(|why| format!("{name}:{line_no}: {why}"))?;
    }
}

fn number(field: &str, what: &str) -> Result<u64, String> {
    field
        .parse()
        .map_err(|_| format!("{what} is not a whole number: {field:?}"))
}

/// The value of attribute `key` in a GTF attribute column.
fn attribute<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|a| {
        let v = a.trim_start().strip_prefix(key)?.strip_prefix(' ')?;
        v.strip_prefix('"')?.strip_suffix('"')
    })
}

/// Parses the gene lines of a GTF. Lines of other features, comments and
/// genes off the map's chromosomes are skipped. Returns the genes sorted by
/// chromosome, start, end and name, or `file:line: reason`.
pub fn parse_gtf(input: impl BufRead, name: &str) -> Result<Vec<Gene>, String> {
    let mut genes = Vec::new();
    each_line(input, name, |line| {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 9 {
            return Err(format!("expected 9 columns, found {}", f.len()));
        }
        if f[2] != "gene" {
            return Ok(());
        }
        let Some(chrom) = chrom_index(f[0]) else {
            return Ok(());
        };
        let start = number(f[3], "start")?;
        let end = number(f[4], "end")?;
        if start == 0 || end < start {
            return Err(format!("not a 1-based span: {start}-{end}"));
        }
        let strand = match f[6] {
            "+" => b'+',
            "-" => b'-',
            "." => b'.',
            s => return Err(format!("strand is not +, - or .: {s:?}")),
        };
        let gene_name = attribute(f[8], "gene_name")
            .or_else(|| attribute(f[8], "gene_id"))
            .ok_or("no gene_name and no gene_id")?;
        genes.push(Gene {
            chrom,
            s: start - 1,
            e: end,
            strand,
            name: gene_name.to_owned(),
            kind: attribute(f[8], "gene_type").unwrap_or("").to_owned(),
        });
        Ok(())
    })?;
    genes.sort_unstable_by(|a, b| {
        (a.chrom, a.s, a.e, &a.name, &a.kind, a.strand)
            .cmp(&(b.chrom, b.s, b.e, &b.name, &b.kind, b.strand))
    });
    Ok(genes)
}

/// Sorts spans and merges each one that overlaps or touches the block
/// before it, so `[10, 20)` and `[20, 30)` become `[10, 30)` and
/// `[10, 20)` and `[21, 30)` stay two blocks.
pub fn merge(spans: &mut Vec<(u64, u64)>) {
    spans.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(spans.len());
    for &(s, e) in spans.iter() {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    *spans = out;
}

/// Parses an exon BED (chrom, 0-based start, end; more columns are
/// ignored) into merged blocks per chromosome, or `file:line: reason`.
pub fn parse_bed(input: impl BufRead, name: &str) -> Result<Vec<Vec<(u64, u64)>>, String> {
    let mut per: Vec<Vec<(u64, u64)>> = vec![Vec::new(); CHROMS.len()];
    each_line(input, name, |line| {
        let mut f = line.split('\t');
        let (Some(c), Some(s), Some(e)) = (f.next(), f.next(), f.next()) else {
            return Err("expected at least 3 columns".to_owned());
        };
        let Some(chrom) = chrom_index(c) else {
            return Ok(());
        };
        let (s, e) = (number(s, "start")?, number(e, "end")?);
        if e <= s {
            return Err(format!("not a BED span: {s}-{e}"));
        }
        per[chrom].push((s, e));
        Ok(())
    })?;
    for spans in &mut per {
        merge(spans);
    }
    Ok(per)
}

/// Reads both inputs.
pub fn load(inputs: &Inputs) -> Result<Genes, BuildError> {
    let read = |path: &Path| -> Result<_, BuildError> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let input =
            og_io::open(path).map_err(|e| BuildError::NoMap(format!("{MISSING}: {name}: {e}")))?;
        Ok((input, name))
    };
    let unreadable = |e: String| BuildError::NoMap(format!("the gene files cannot be built: {e}"));
    let (gtf, gtf_name) = read(&inputs.gtf)?;
    let genes = parse_gtf(gtf, &gtf_name).map_err(unreadable)?;
    let (bed, bed_name) = read(&inputs.bed)?;
    let exons = parse_bed(bed, &bed_name).map_err(unreadable)?;
    Ok(Genes {
        release: inputs.release.clone(),
        genes,
        exons,
    })
}

fn ints(xs: impl IntoIterator<Item = u64>) -> Value {
    Value::Array(xs.into_iter().map(|x| Value::Int(x.into())).collect())
}

fn strs<'a>(xs: impl IntoIterator<Item = &'a str>) -> Value {
    Value::Array(xs.into_iter().map(|x| Value::Str(x.to_owned())).collect())
}

/// One chromosome's file as a JSON value.
pub fn chrom_value(g: &Genes, chrom: usize) -> Value {
    let lo = g.genes.partition_point(|x| x.chrom < chrom);
    let hi = g.genes.partition_point(|x| x.chrom <= chrom);
    let genes = &g.genes[lo..hi];
    let blocks = &g.exons[chrom];
    obj([
        ("v", Value::Int(VERSION.into())),
        ("gencode", Value::Str(g.release.clone())),
        ("chrom", Value::Str(CHROMS[chrom].to_owned())),
        (
            "g",
            obj([
                ("ds", ints(gaps(genes.iter().map(|x| x.s)))),
                ("l", ints(genes.iter().map(|x| x.e - x.s))),
                (
                    "strand",
                    strs(genes.iter().map(|x| match x.strand {
                        b'+' => "+",
                        b'-' => "-",
                        _ => ".",
                    })),
                ),
                ("name", strs(genes.iter().map(|x| x.name.as_str()))),
                ("type", strs(genes.iter().map(|x| x.kind.as_str()))),
            ]),
        ),
        (
            "x",
            obj([
                ("ds", ints(gaps(blocks.iter().map(|b| b.0)))),
                ("l", ints(blocks.iter().map(|b| b.1 - b.0))),
            ]),
        ),
    ])
}

/// The names file as a JSON value.
pub fn names_value(g: &Genes) -> Value {
    let mut by_name: Vec<&Gene> = g.genes.iter().collect();
    by_name.sort_by(|a, b| (&a.name, a.chrom, a.s, a.e).cmp(&(&b.name, b.chrom, b.s, b.e)));
    obj([
        ("v", Value::Int(VERSION.into())),
        ("gencode", Value::Str(g.release.clone())),
        ("name", strs(by_name.iter().map(|x| x.name.as_str()))),
        ("chrom", strs(by_name.iter().map(|x| CHROMS[x.chrom]))),
        ("s", ints(by_name.iter().map(|x| x.s))),
        ("e", ints(by_name.iter().map(|x| x.e))),
    ])
}

/// `value` as compact JSON in one gzip member with mtime 0 and no name.
fn gzip_json(value: &Value) -> io::Result<Vec<u8>> {
    let opts = Dump {
        separators: Some((",", ":")),
        ..Dump::default()
    };
    let text =
        json::dumps(value, &opts).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut gz = GzBuilder::new()
        .mtime(0)
        .write(Vec::with_capacity(text.len() / 4), Compression::best());
    gz.write_all(text.as_bytes())?;
    gz.finish()
}

/// Every file of the set, as (file name, bytes): the 24 chromosomes, then
/// the names. One worker per core takes the next file until none is left,
/// so at most that many files are in memory as JSON at once.
pub fn write_all(g: &Genes) -> io::Result<Vec<(String, Vec<u8>)>> {
    let wants: Vec<Want> = (0..CHROMS.len())
        .map(Want::Chrom)
        .chain([Want::Names])
        .collect();
    let next = AtomicUsize::new(0);
    let workers = std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .min(wants.len());
    let mut done: Vec<(usize, io::Result<Vec<u8>>)> = std::thread::scope(|scope| {
        let jobs: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(&w) = wants.get(i) else {
                            return out;
                        };
                        let value = match w {
                            Want::Chrom(c) => chrom_value(g, c),
                            Want::Names => names_value(g),
                        };
                        out.push((i, gzip_json(&value)));
                    }
                })
            })
            .collect();
        jobs.into_iter()
            .flat_map(|j| j.join().expect("a gene file worker does not panic"))
            .collect()
    });
    done.sort_unstable_by_key(|(i, _)| *i);
    done.into_iter()
        .map(|(i, bytes)| bytes.map(|b| (wants[i].file_name(), b)))
        .collect()
}

/// The cache folder of one input stamp.
pub fn cache_dir(cache: &Path, stamp: &str) -> PathBuf {
    cache.join("genes").join(stamp)
}

/// Moves a finished set into place, then removes older sets. A set that
/// another writer moved in first is kept and ours dropped.
fn store(dir: &Path, files: &[(String, Vec<u8>)]) -> io::Result<()> {
    let parent = dir.parent().expect("a cache folder has a parent");
    let stamp = dir.file_name().expect("a cache folder has a name");
    fs::create_dir_all(parent)?;
    let tmp = parent.join(format!(
        ".{}.tmp-{}",
        stamp.to_string_lossy(),
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&tmp);
    fs::create_dir(&tmp)?;
    let written = files
        .iter()
        .try_for_each(|(name, bytes)| fs::write(tmp.join(name), bytes));
    if let Err(e) = written.and_then(|()| fs::rename(&tmp, dir)) {
        let _ = fs::remove_dir_all(&tmp);
        if !dir.is_dir() {
            return Err(e);
        }
    }
    for entry in fs::read_dir(parent)? {
        let entry = entry?;
        let name = entry.file_name();
        // Another writer's set in progress starts with a dot: never touch it.
        if name != stamp && !name.to_string_lossy().starts_with('.') && entry.path().is_dir() {
            fs::remove_dir_all(entry.path())?;
        }
    }
    Ok(())
}

/// The bytes of the file `want` asks for: from the cache when the inputs
/// have not changed, else built in full, stored and returned. A cache that
/// cannot be written is reported on stderr and the file is still returned.
pub fn genes_cached(
    root: &Path,
    cache: &Path,
    made_by: &str,
    want: Want,
) -> Result<Vec<u8>, BuildError> {
    let inputs = locate(root)?;
    let stamp = crate::stamp::stamp(
        VERSION,
        made_by,
        b"genes",
        &[inputs.gtf.clone(), inputs.bed.clone()],
    );
    let dir = cache_dir(cache, &stamp);
    let file = want.file_name();
    if let Ok(bytes) = fs::read(dir.join(&file)) {
        return Ok(bytes);
    }
    let genes = load(&inputs)?;
    let files = write_all(&genes).map_err(|e| BuildError::Internal(e.to_string()))?;
    if let Err(e) = store(&dir, &files) {
        eprintln!("og map genes: the cache was not written: {e}");
    }
    files
        .into_iter()
        .find(|(name, _)| *name == file)
        .map(|(_, bytes)| bytes)
        .ok_or_else(|| BuildError::Internal(format!("no file {file} in the set")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gtf(lines: &str) -> Result<Vec<Gene>, String> {
        parse_gtf(lines.as_bytes(), "genes.gtf")
    }

    #[test]
    fn gtf_genes_become_zero_based_spans() {
        let text = "##description: test\n\
chr2\tHAVANA\tgene\t100\t200\t.\t-\t.\tgene_id \"ENSG2.1\"; gene_type \"lncRNA\"; gene_name \"B\"; level 2;\n\
chr2\tHAVANA\ttranscript\t100\t200\t.\t-\t.\tgene_id \"ENSG2.1\";\n\
chr1\tENSEMBL\tgene\t1\t1\t.\t+\t.\tgene_id \"ENSG1.1\"; gene_type \"protein_coding\"; gene_name \"A\";\n\
chrM\tENSEMBL\tgene\t648\t1601\t.\t+\t.\tgene_id \"ENSG3.1\"; gene_name \"MT-RNR1\";\n\
chr1_KI270706v1_random\tENSEMBL\tgene\t5\t9\t.\t+\t.\tgene_id \"ENSG4.1\";\n\
chrX\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id \"ENSG5.1\";\n";
        let g = gtf(text).unwrap();
        let rows: Vec<_> = g
            .iter()
            .map(|x| {
                (
                    CHROMS[x.chrom],
                    x.s,
                    x.e,
                    x.strand,
                    x.name.as_str(),
                    x.kind.as_str(),
                )
            })
            .collect();
        assert_eq!(
            rows,
            [
                ("chr1", 0, 1, b'+', "A", "protein_coding"),
                ("chr2", 99, 200, b'-', "B", "lncRNA"),
                ("chrX", 9, 20, b'+', "ENSG5.1", ""),
            ]
        );
    }

    #[test]
    fn gtf_errors_name_the_line() {
        let bad = "chr1\tH\tgene\t0\t5\t.\t+\t.\tgene_id \"X\";\n";
        assert_eq!(
            gtf(bad).unwrap_err(),
            "genes.gtf:1: not a 1-based span: 0-5"
        );
        let short = "\nchr1\tH\tgene\t1\t5\n";
        assert_eq!(
            gtf(short).unwrap_err(),
            "genes.gtf:2: expected 9 columns, found 5"
        );
        let strand = "chr1\tH\tgene\t1\t5\t.\t?\t.\tgene_id \"X\";\n";
        assert!(gtf(strand).unwrap_err().contains("strand"));
        let anon = "chr1\tH\tgene\t1\t5\t.\t+\t.\tlevel 2;\n";
        assert_eq!(
            gtf(anon).unwrap_err(),
            "genes.gtf:1: no gene_name and no gene_id"
        );
    }

    #[test]
    fn exons_merge_when_they_overlap_or_touch() {
        let mut v = vec![(30, 40), (10, 20), (20, 25), (12, 15), (26, 28), (39, 45)];
        merge(&mut v);
        assert_eq!(v, [(10, 25), (26, 28), (30, 45)]);
        let mut gap = vec![(10, 20), (21, 30)];
        merge(&mut gap);
        assert_eq!(gap, [(10, 20), (21, 30)]);
        let mut empty: Vec<(u64, u64)> = Vec::new();
        merge(&mut empty);
        assert!(empty.is_empty());
    }

    #[test]
    fn bed_rows_merge_per_chromosome() {
        let text = "chr1\t10\t20\nchr1\t15\t30\nchr2\t30\t40\nchrM\t1\t9\nchr1\t30\t31\n";
        let per = parse_bed(text.as_bytes(), "exons.bed").unwrap();
        assert_eq!(per[0], [(10, 31)]);
        assert_eq!(per[1], [(30, 40)]);
        assert!(per[2..].iter().all(Vec::is_empty));
        assert_eq!(
            parse_bed("chr1\t5\t5\n".as_bytes(), "exons.bed").unwrap_err(),
            "exons.bed:1: not a BED span: 5-5"
        );
        assert_eq!(
            parse_bed("chr1\t5\n".as_bytes(), "exons.bed").unwrap_err(),
            "exons.bed:1: expected at least 3 columns"
        );
    }

    #[test]
    fn release_from_the_file_name() {
        assert_eq!(release_of("gencode.v50.basic.genes.gtf"), Some(50));
        assert_eq!(release_of("gencode.v7.basic.genes.gtf"), Some(7));
        assert_eq!(release_of("gencode.v.basic.genes.gtf"), None);
        assert_eq!(release_of("gencode.v50.basic.exons.bed"), None);
        assert_eq!(release_of("gencode.v5x.basic.genes.gtf"), None);
    }

    #[test]
    fn want_names_its_file() {
        assert_eq!(
            Want::chrom("chr7").unwrap().file_name(),
            "genes.chr7.json.gz"
        );
        assert_eq!(Want::Names.file_name(), "genes.names.json.gz");
        assert!(matches!(Want::chrom("chrM"), Err(BuildError::Usage(_))));
    }
}
