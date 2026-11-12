//! `og mito annotate|deletions|copies`: step 20's chrM tables.
//!
//! - `annotate`: each ALT of the filtered chrM VCF with its ClinVar record
//!   (CLNSIG, CLNREVSTAT, CLNDN) and gnomAD v3.1's chrM sites (AF_hom,
//!   AF_het, max_hl), matched on the allele trimmed to its shortest form
//!   ([`min_allele`]), so Mutect2's, ClinVar's and gnomAD's spellings of one
//!   indel meet. `<id>_chrM_annotated.tsv`.
//! - `deletions`: eKLIPse's deletions CSV (Goudenege 2019) as a table; the
//!   deleted length wraps the circular origin. `<id>_chrM_deletions.tsv`.
//! - `copies`: the mtDNA copy number, 2 x (chrM reads per base) / (autosomal
//!   reads per base) from `samtools idxstats`, not assessed under 20 reads
//!   per kb of autosome. `<id>_mtdna_copies.tsv`.

use crate::args::{self, Kind, Opt};
use crate::{
    basename, finish, float, int, is_file, read_text, repr, text_lines, trunc, write_file,
};
use og_io::tabix;
use pycompat::fmt::{repr as float_repr, round};
use std::collections::HashMap;
use std::io::Write;

/// The table columns of `annotate`.
pub const ANNOT_COLS: [&str; 13] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "filter",
    "af",
    "clinvar_id",
    "clinvar_significance",
    "clinvar_review",
    "clinvar_condition",
    "gnomad_af_hom",
    "gnomad_af_het",
    "gnomad_max_hl",
];
/// The table columns of `deletions`.
pub const DEL_COLS: [&str; 11] = [
    "bkp5",
    "bkp3",
    "length",
    "level_pct",
    "level_fwd_pct",
    "level_rev_pct",
    "blast_5",
    "blast_3",
    "depth_5",
    "depth_3",
    "repeat",
];
/// The rows of `copies`.
pub const CN_COLS: [&str; 8] = [
    "status",
    "reason",
    "chrm_reads",
    "chrm_length",
    "autosomal_reads",
    "autosomal_length",
    "depth_ratio",
    "mtdna_copies",
];
/// The rCRS, GRCh38's chrM.
pub const MITO_LENGTH: i128 = 16569;
/// Under this many reads per kb of autosome the ratio means nothing.
pub const CN_MIN_READS_PER_KB: f64 = 20.0;
/// The names a BAM gives chrM.
const CHRM_NAMES: [&str; 4] = ["chrM", "MT", "chrMT", "M"];

/// One allele as a key: (pos, ref, alt).
pub type Allele = (i128, String, String);

/// `(pos, ref, alt)` without the bases ref and alt share at the end and then
/// at the start, one base always kept: Mutect2's `A>AC` of a two-ALT record
/// `AC>A,ACC` becomes the same key as ClinVar's and gnomAD's.
pub fn min_allele(mut pos: i128, r: &str, a: &str) -> Allele {
    let (mut r, mut a): (Vec<char>, Vec<char>) = (r.chars().collect(), a.chars().collect());
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    let mut k = 0;
    while r.len() - k > 1 && a.len() - k > 1 && r[k] == a[k] {
        k += 1;
        pos += 1;
    }
    (pos, r[k..].iter().collect(), a[k..].iter().collect())
}

/// `{k: v}` of a VCF INFO field; a key without `=` maps to `""`.
fn info_map(field: &str) -> HashMap<&str, &str> {
    field
        .split(';')
        .map(|kv| kv.split_once('=').unwrap_or((kv, "")))
        .collect()
}

/// The split data lines of a VCF, plain or gzip.
fn vcf_records(path: &str) -> Result<Vec<Vec<String>>, String> {
    let lines = text_lines(path).map_err(|e| format!("{}: {e}", basename(path)))?;
    let mut out = Vec::new();
    for line in lines {
        let line = line.map_err(|e| format!("{}: {e}", basename(path)))?;
        if !line.starts_with('#') {
            out.push(line.split('\t').map(str::to_owned).collect());
        }
    }
    Ok(out)
}

/// The `##fileDate=` of a VCF's header.
fn vcf_file_date(path: &str) -> Result<Option<String>, String> {
    for line in text_lines(path)? {
        let line = line.map_err(|e| format!("{path}: {e}"))?;
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("##fileDate=") {
            let v = pycompat::text::strip(&line, None);
            return Ok(Some(v.split_once('=').map_or("", |x| x.1).to_owned()));
        }
    }
    Ok(None)
}

/// ClinVar's columns for an allele: id, CLNSIG, CLNREVSTAT, CLNDN.
type ClinVar = [String; 4];
/// gnomAD's columns for an allele: AF_hom, AF_het, max_hl.
type Gnomad = [String; 3];

/// The ClinVar records of chrM at the positions `wanted`, each read and its
/// alleles keyed in the order the Python's lookups visited them: for each
/// wanted position in order, the records starting a base before, at and a
/// base after it (ClinVar's record can start a base away from Mutect2's for
/// the same indel), the last visit winning.
fn clinvar_alleles(path: &str, wanted: &[i128]) -> Result<HashMap<Allele, ClinVar>, String> {
    let mut out = HashMap::new();
    let mut reader = tabix::Reader::open(path).map_err(|e| format!("{}: {e}", basename(path)))?;
    let Some(chrom) = ["chrM", "MT"]
        .into_iter()
        .find(|c| reader.index().names().iter().any(|n| n == c))
    else {
        return Ok(out);
    };
    // chrM's records once, by their start.
    let mut at: HashMap<i128, Vec<Vec<String>>> = HashMap::new();
    let last = wanted.iter().copied().max().unwrap_or(0) + 1;
    let end = u64::try_from(last.max(1)).unwrap_or(u64::MAX).min(1 << 29);
    for line in reader.query(chrom, 0, end) {
        let line = line.map_err(|e| format!("{}: {e}", basename(path)))?;
        let line = String::from_utf8_lossy(&line).into_owned();
        let f: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if f.len() < 3 || f[1].is_empty() || !f[1].bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }
        let pos = int(&f[1])?;
        at.entry(pos).or_default().push(f);
    }
    for &pos in wanted {
        for q in [pos - 1, pos, pos + 1] {
            if q < 1 {
                continue;
            }
            for c in at.get(&q).into_iter().flatten() {
                let i = c.get(7).map(|x| info_map(x)).unwrap_or_default();
                let g = |k: &str| i.get(k).copied().unwrap_or(".").to_owned();
                let alts = c.get(4).ok_or_else(|| {
                    format!("{}: a record with {} fields", basename(path), c.len())
                })?;
                for alt in alts.split(',') {
                    out.insert(
                        min_allele(q, &c[3], alt),
                        [c[2].clone(), g("CLNSIG"), g("CLNREVSTAT"), g("CLNDN")],
                    );
                }
            }
        }
    }
    Ok(out)
}

/// The annotation table: its two `#` header lines and its rows.
pub fn annotate(
    vcf: &str,
    clinvar: Option<&str>,
    gnomad: Option<&str>,
) -> Result<(Vec<String>, Vec<[String; 13]>), String> {
    let mut head = Vec::new();
    let clinvar = match clinvar {
        Some(c) if is_file(c) && is_file(format!("{c}.tbi")) => {
            let date = vcf_file_date(c)?.unwrap_or_else(|| "fileDate unknown".into());
            head.push(format!("# clinvar\t{}\t{date}", basename(c)));
            Some(c)
        }
        Some(c) => {
            head.push(format!("# clinvar\tnot annotated\tno {c} with its .tbi"));
            None
        }
        None => {
            head.push("# clinvar\tnot annotated\tno ClinVar file given".into());
            None
        }
    };
    let mut gn: HashMap<Allele, Gnomad> = HashMap::new();
    match gnomad {
        Some(g) if is_file(g) => {
            head.push(format!("# gnomad\t{}\t-", basename(g)));
            for r in vcf_records(g)? {
                if r.len() < 8 {
                    continue;
                }
                let i = info_map(&r[7]);
                let v = |k: &str| i.get(k).copied().unwrap_or(".").to_owned();
                let pos = int(&r[1])?;
                for alt in r[4].split(',') {
                    gn.insert(
                        min_allele(pos, &r[3], alt),
                        [v("AF_hom"), v("AF_het"), v("max_hl")],
                    );
                }
            }
        }
        Some(g) => head.push(format!("# gnomad\tnot annotated\tno {g}")),
        None => head.push("# gnomad\tnot annotated\tno gnomAD chrM file given".into()),
    }
    let recs = vcf_records(vcf)?;
    let mut cv = HashMap::new();
    if let Some(c) = clinvar {
        let mut wanted = Vec::new();
        for r in recs.iter().filter(|r| r.len() >= 5) {
            let pos = int(&r[1])?;
            wanted.push(pos);
            for a in r[4].split(',') {
                wanted.push(min_allele(pos, &r[3], a).0);
            }
        }
        wanted.sort_unstable();
        wanted.dedup();
        cv = clinvar_alleles(c, &wanted)?;
    }
    let dot = || ".".to_owned();
    let mut rows = Vec::new();
    for r in recs.iter().filter(|r| r.len() >= 8) {
        let (fmt, val): (Vec<&str>, Vec<&str>) = if r.len() > 9 {
            (r[8].split(':').collect(), r[9].split(':').collect())
        } else {
            (Vec::new(), Vec::new())
        };
        let afs: Vec<&str> = match fmt.iter().position(|f| *f == "AF") {
            Some(k) if k < val.len() => val[k].split(',').collect(),
            _ => Vec::new(),
        };
        let pos = int(&r[1])?;
        for (k, alt) in r[4].split(',').enumerate() {
            let key = min_allele(pos, &r[3], alt);
            let [cid, sig, rev, dn] = cv
                .get(&key)
                .cloned()
                .unwrap_or_else(|| [dot(), dot(), dot(), dot()]);
            let [hom, het, hl] = gn
                .get(&key)
                .cloned()
                .unwrap_or_else(|| [dot(), dot(), dot()]);
            rows.push([
                r[0].clone(),
                key.0.to_string(),
                key.1,
                key.2,
                r[6].clone(),
                afs.get(k).map_or_else(dot, |s| (*s).to_owned()),
                cid,
                sig,
                rev,
                dn,
                hom,
                het,
                hl,
            ]);
        }
    }
    Ok((head, rows))
}

/// `float(x.replace(",", "."))` of eKLIPse's numbers, finite only.
fn eklipse_num(x: &str) -> Result<f64, String> {
    let v = float(&x.replace(',', "."))?;
    if !v.is_finite() {
        return Err(format!("{} is not a number", repr(x)));
    }
    Ok(v)
}

/// One deletion of eKLIPse's table.
#[derive(Debug, Clone, PartialEq)]
pub struct Deletion {
    pub bkp5: i128,
    pub bkp3: i128,
    pub length: i128,
    pub level_pct: f64,
    pub level_fwd_pct: f64,
    pub level_rev_pct: f64,
    pub blast_5: i128,
    pub blast_3: i128,
    pub depth_5: i128,
    pub depth_3: i128,
    pub repeat: String,
}

impl Deletion {
    fn cells(&self) -> [String; 11] {
        [
            self.bkp5.to_string(),
            self.bkp3.to_string(),
            self.length.to_string(),
            float_repr(self.level_pct),
            float_repr(self.level_fwd_pct),
            float_repr(self.level_rev_pct),
            self.blast_5.to_string(),
            self.blast_3.to_string(),
            self.depth_5.to_string(),
            self.depth_3.to_string(),
            self.repeat.clone(),
        ]
    }
}

/// The repr of a list of strings, as Python prints `rows[0][:3]`.
fn list_repr(v: &[String]) -> String {
    let items: Vec<String> = v.iter().map(|s| repr(s)).collect();
    format!("[{}]", items.join(", "))
}

/// The rows of eKLIPse's `eKLIPse_deletions.csv`: `;` between fields, each
/// quoted, a comma for the decimal point.
pub fn deletions(path: &str) -> Result<Vec<Deletion>, String> {
    let bytes = std::fs::read(path).map_err(|e| crate::os_error(&e, path))?;
    let text = String::from_utf8_lossy(&bytes);
    // newline="": lines end at \n, \r\n or a lone \r, and the csv reader
    // sees the endings as written.
    let mut lines = Vec::new();
    let (b, mut start) = (text.as_bytes(), 0);
    for k in 0..b.len() {
        let ends = b[k] == b'\n' || (b[k] == b'\r' && b.get(k + 1) != Some(&b'\n'));
        if ends {
            lines.push(Ok::<_, std::io::Error>(&text[start..=k]));
            start = k + 1;
        }
    }
    if start < b.len() {
        lines.push(Ok(&text[start..]));
    }
    let lines = lines.into_iter();
    let mut rows = Vec::new();
    for r in pycompat::csv::Reader::new(lines, ';', pycompat::csv::Quoting::Minimal) {
        let r = r.map_err(|e| e.to_string())?;
        if !r.is_empty() {
            rows.push(r);
        }
    }
    let want = ["Title", "5' breakpoint", "3' breakpoint"];
    if rows
        .first()
        .is_none_or(|h| h.iter().take(3).map(String::as_str).ne(want))
    {
        let got = rows
            .first()
            .map_or_else(|| "missing".to_owned(), |h| list_repr(&h[..h.len().min(3)]));
        return Err(format!(
            "{path}: not eKLIPse's deletions table (header {got})"
        ));
    }
    let mut out = Vec::new();
    for r in &rows[1..] {
        if r.len() < 11 {
            return Err(format!(
                "{path}: a row with {} fields, want 11: {}",
                r.len(),
                r.join(";")
            ));
        }
        let level = |k: usize| -> Result<f64, String> {
            round(eklipse_num(&r[k])?, 3).map_err(|e| e.to_string())
        };
        let (b5, b3) = (int(&r[1])?, int(&r[2])?);
        out.push(Deletion {
            bkp5: b5,
            bkp3: b3,
            length: if b3 > b5 {
                b3 - b5 - 1
            } else {
                MITO_LENGTH - b5 + b3 - 1
            },
            level_pct: level(3)?,
            level_fwd_pct: level(4)?,
            level_rev_pct: level(5)?,
            blast_5: int(&r[6])?,
            blast_3: int(&r[7])?,
            depth_5: trunc(eklipse_num(&r[8])?)?,
            depth_3: trunc(eklipse_num(&r[9])?)?,
            repeat: if r[10].is_empty() {
                "none".into()
            } else {
                r[10].clone()
            },
        });
    }
    Ok(out)
}

/// Why the copy number was not assessed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CnReason {
    /// No chrM in the index, or no reads or no length on it.
    NoChrmReads,
    /// Under [`CN_MIN_READS_PER_KB`] on the autosomes.
    LowAutosomalDepth,
}

impl CnReason {
    fn word(self) -> &'static str {
        match self {
            CnReason::NoChrmReads => "no_chrM_reads",
            CnReason::LowAutosomalDepth => "low_autosomal_depth",
        }
    }
}

/// The copy number block.
#[derive(Debug, Clone, PartialEq)]
pub struct Copies {
    pub chrm: Option<(i128, i128)>,
    pub autosomal_reads: i128,
    pub autosomal_length: i128,
    /// `(depth ratio, copies)` when called, else why not.
    pub call: Result<(f64, f64), CnReason>,
}

impl Copies {
    fn rows(&self) -> Vec<(&'static str, String)> {
        let dot = || ".".to_owned();
        let (ratio, copies) = match self.call {
            Ok((r, c)) => (float_repr(r), float_repr(c)),
            Err(_) => (dot(), dot()),
        };
        vec![
            (
                "status",
                if self.call.is_ok() {
                    "called"
                } else {
                    "not assessed"
                }
                .to_owned(),
            ),
            (
                "reason",
                self.call.err().map_or_else(dot, |r| r.word().to_owned()),
            ),
            (
                "chrm_reads",
                self.chrm.map_or_else(dot, |c| c.0.to_string()),
            ),
            (
                "chrm_length",
                self.chrm.map_or_else(dot, |c| c.1.to_string()),
            ),
            ("autosomal_reads", self.autosomal_reads.to_string()),
            ("autosomal_length", self.autosomal_length.to_string()),
            ("depth_ratio", ratio),
            ("mtdna_copies", copies),
        ]
    }
}

/// True for chr1-chr22 and 1-22.
fn autosome(name: &str) -> bool {
    let n = name.strip_prefix("chr").unwrap_or(name);
    !n.is_empty()
        && n.len() <= 2
        && !n.starts_with('0')
        && n.bytes().all(|b| b.is_ascii_digit())
        && n.parse::<u8>().is_ok_and(|v| (1..=22).contains(&v))
}

/// `str.splitlines()`: the line boundaries Python knows.
fn splitlines(s: &str) -> impl Iterator<Item = &str> {
    let mut parts: Vec<&str> = s
        .split([
            '\n', '\r', '\x0b', '\x0c', '\x1c', '\x1d', '\x1e', '\u{85}', '\u{2028}', '\u{2029}',
        ])
        .collect();
    if parts.last() == Some(&"") {
        parts.pop();
    }
    parts.into_iter()
}

/// The copy number from `samtools idxstats` text: chrM reads per base over
/// autosomal reads per base, times 2.
pub fn copies(idxstats: &str) -> Result<Copies, String> {
    let mut chrm = None;
    let (mut auto_reads, mut auto_len) = (0i128, 0i128);
    for line in splitlines(idxstats) {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 3 || c[0] == "*" {
            continue;
        }
        let (length, mapped) = (int(c[1])?, int(c[2])?);
        if CHRM_NAMES.contains(&c[0]) {
            chrm = Some((mapped, length));
        } else if autosome(c[0]) {
            auto_reads += mapped;
            auto_len += length;
        }
    }
    let call = match chrm {
        Some((r, l)) if r != 0 && l != 0 => {
            if auto_len == 0 || (auto_reads as f64 / auto_len as f64) * 1000.0 < CN_MIN_READS_PER_KB
            {
                Err(CnReason::LowAutosomalDepth)
            } else {
                let ratio = (r as f64 / l as f64) / (auto_reads as f64 / auto_len as f64);
                let r3 = round(ratio, 3).map_err(|e| e.to_string())?;
                let c1 = round(2.0 * ratio, 1).map_err(|e| e.to_string())?;
                Ok((r3, c1))
            }
        }
        _ => Err(CnReason::NoChrmReads),
    };
    Ok(Copies {
        chrm,
        autosomal_reads: auto_reads,
        autosomal_length: auto_len,
        call,
    })
}

const ANNOTATE_USAGE: &str =
    "usage: og mito annotate --vcf VCF [--clinvar VCF] [--gnomad VCF] --out F\n";
const DELETIONS_USAGE: &str = "usage: og mito deletions --eklipse CSV --out F\n";
const COPIES_USAGE: &str = "usage: og mito copies --idxstats F --out F\n";

/// `og mito VERB ARGS`: the exit code.
pub fn cli(verb: Option<&str>, argv: &[String], out: &mut impl Write, err: &mut impl Write) -> u8 {
    let (usage, opts): (&str, &[Opt]) = match verb {
        Some("annotate") => (
            ANNOTATE_USAGE,
            &[
                Opt("--vcf", Kind::Value, true),
                Opt("--clinvar", Kind::Value, false),
                Opt("--gnomad", Kind::Value, false),
                Opt("--out", Kind::Value, true),
            ],
        ),
        Some("deletions") => (
            DELETIONS_USAGE,
            &[
                Opt("--eklipse", Kind::Value, true),
                Opt("--out", Kind::Value, true),
            ],
        ),
        Some("copies") => (
            COPIES_USAGE,
            &[
                Opt("--idxstats", Kind::Value, true),
                Opt("--out", Kind::Value, true),
            ],
        ),
        _ => {
            let _ = writeln!(err, "usage: og mito annotate|deletions|copies [args]");
            let _ = writeln!(
                err,
                "og mito: error: og mito takes annotate, deletions or copies"
            );
            return 2;
        }
    };
    let verb = verb.unwrap_or_default();
    let a = match args::parse(argv, opts) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{usage}");
            let _ = writeln!(err, "og mito {verb}: error: {e}");
            return 2;
        }
    };
    let dest = a.text("--out");
    let mut say = |s: String| writeln!(out, "{s}").map_err(|e| e.to_string());
    let mut run = || -> Result<(), String> {
        match verb {
            "annotate" => {
                let some = |k: &str| a.get(k).filter(|v| !v.is_empty());
                let (head, rows) = annotate(a.text("--vcf"), some("--clinvar"), some("--gnomad"))?;
                let mut body = String::new();
                for h in &head {
                    body.push_str(h);
                    body.push('\n');
                }
                body.push_str(&ANNOT_COLS.join("\t"));
                body.push('\n');
                for r in &rows {
                    body.push_str(&r.join("\t"));
                    body.push('\n');
                }
                write_file(dest, body.as_bytes())?;
                for h in &head {
                    let c: Vec<&str> = h[2..].split('\t').collect();
                    let note = if c[2] == "-" {
                        String::new()
                    } else {
                        format!(" ({})", c[2])
                    };
                    say(format!("  {}: {}{note}", c[0], c[1]))?;
                }
                let n_cv = rows.iter().filter(|r| r[7] != ".").count();
                let n_gn = rows.iter().filter(|r| r[10] != ".").count();
                say(format!(
                    "  {} chrM alleles: {n_cv} with a ClinVar record, {n_gn} in gnomAD",
                    rows.len()
                ))
            }
            "deletions" => {
                let rows = deletions(a.text("--eklipse"))?;
                let mut body = DEL_COLS.join("\t");
                body.push('\n');
                for r in &rows {
                    body.push_str(&r.cells().join("\t"));
                    body.push('\n');
                }
                write_file(dest, body.as_bytes())?;
                say(format!("  {} large deletions (eKLIPse)", rows.len()))
            }
            _ => {
                let b = copies(&read_text(a.text("--idxstats"))?)?;
                let mut body = String::from("metric\tvalue\n");
                for (k, v) in b.rows() {
                    body.push_str(&format!("{k}\t{v}\n"));
                }
                write_file(dest, body.as_bytes())?;
                let shown = match b.call {
                    Ok((_, c)) => float_repr(c),
                    Err(r) => format!("not assessed ({})", r.word()),
                };
                say(format!("  mtDNA copy number: {shown}"))
            }
        }
    };
    finish(run(), err)
}

#[cfg(test)]
mod tests;
