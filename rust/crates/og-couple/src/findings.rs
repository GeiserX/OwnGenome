//! Each partner's carrier findings, one reader per step output: step 06's
//! ClinVar hits, step 19's Delly deletions, step 09's repeats, step 35's
//! SMN1 copy number and CYP21A2 table.

use crate::files::{each_line, each_record, info_get, int, is_digits};
use crate::list::{CarrierList, GeneIndex};
use og_io::Utf8;
use og_pedigree::{Die, Par};
use og_steps::sv::{DEPTH_MIN_SIZE, Rdcn, Sex, SvKind, depth_contradicts};
use og_summary::collect::clinvar::{Floor, Hit};
use std::collections::HashMap;

/// A structural variant is 50 bp or longer. A shorter Delly deletion is
/// left to the small-variant caller and counts only through ClinVar.
pub const SV_MIN_SIZE: i64 = 50;
/// A benign loss record of 1% or more (AnnotSV's own `-benignAF` default)
/// is a common deletion.
pub const COMMON_AF: f64 = 0.01;
/// A Delly deletion matches a common one when each covers at least half
/// of the other.
pub const COMMON_OVERLAP_PERCENT: i64 = 50;
/// A Delly deletion longer than this needs read depth to show the loss
/// (`og_steps::sv::depth_contradicts`).
pub const DEPTH_MIN: i64 = DEPTH_MIN_SIZE;
/// FMR1: 55 to 200 CGG is a premutation, more than 200 a full mutation.
const FMR1_FULL: i64 = 201;
/// The repeat loci of the list: the smallest repeat count counted as a
/// carrier finding, and the repeat unit.
const REPEATS: [(&str, i64, &str); 2] = [("FMR1", 55, "CGG"), ("FXN", 66, "GAA")];

/// One carrier finding of one person in one gene.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Finding {
    pub gene: String,
    pub text: String,
    /// On how many of the person's copies: 1 or 2.
    pub copies: u8,
    /// The name of a low-function allele of the list's severity column: it
    /// counts only beside a full loss of function in the partner.
    pub hypomorph: Option<String>,
}

impl Finding {
    pub fn new(gene: &str, text: String, copies: u8) -> Finding {
        Finding {
            gene: gene.to_owned(),
            text,
            copies,
            hypomorph: None,
        }
    }
}

/// How many copies a call carries, as the reports word it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    Het,
    Hom,
    Hemizygous,
}

impl Zygosity {
    pub fn word(self) -> &'static str {
        match self {
            Zygosity::Het => "het",
            Zygosity::Hom => "hom",
            Zygosity::Hemizygous => "hemizygous",
        }
    }

    /// The copies a finding of this zygosity is on.
    pub fn copies(self) -> u8 {
        if self == Zygosity::Het { 1 } else { 2 }
    }
}

fn sex_chrom(chrom: &str) -> bool {
    matches!(chrom, "chrX" | "chrY")
}

/// The zygosity of a GT, or `None` when it carries no alt allele. A `.`
/// allele still counts in the ploidy: step 06 splits a GT 1/2 record into
/// `./1` halves, and each half is one copy, never hom.
pub fn zygosity(gt: &str, chrom: &str, sex: Option<&str>, in_par: bool) -> Option<Zygosity> {
    let alleles: Vec<&str> = gt.split(['/', '|']).filter(|a| !a.is_empty()).collect();
    let alt = alleles.iter().filter(|a| !matches!(**a, "0" | ".")).count();
    if alt == 0 {
        return None;
    }
    if sex_chrom(chrom) && sex == Some("male") && !in_par {
        return Some(Zygosity::Hemizygous);
    }
    if alleles.len() == 1 {
        return Some(if sex_chrom(chrom) {
            Zygosity::Hemizygous
        } else {
            Zygosity::Hom
        });
    }
    Some(if alt == alleles.len() {
        Zygosity::Hom
    } else {
        Zygosity::Het
    })
}

/// One alt allele and one missing one: `1/.`, `./1`, `1|.`.
fn half_call(gt: &str) -> bool {
    let b = gt.as_bytes();
    if b.len() != 3 || !matches!(b[1], b'/' | b'|') {
        return false;
    }
    let alt = |c: u8| (b'1'..=b'9').contains(&c);
    (alt(b[0]) && b[2] == b'.') || (b[0] == b'.' && alt(b[2]))
}

/// `1 star`, `2 stars`.
fn stars_text(n: u8) -> String {
    format!("{n} star{}", if n == 1 { "" } else { "s" })
}

/// Step 06's hits of one person.
pub struct Clinvar {
    pub findings: Vec<Finding>,
    /// Every hit, in a listed gene or not.
    pub total: usize,
    /// The hits the list marks as hypomorph or pseudodeficiency: no carrier
    /// finding on their own.
    pub aside: Vec<String>,
}

/// Step 06's hits in genes of the list. A hit the list marks as a
/// hypomorph or pseudodeficiency allele is set aside as a text; a
/// hypomorph is also kept as a finding marked hypomorph, which beside a
/// full loss of function in the partner's same gene is the child's second
/// hit (TAR syndrome is a null allele plus an RBM8A low-expression allele).
pub fn clinvar(
    path: &str,
    list: &CarrierList,
    sex: Option<&str>,
    par: &Par,
) -> Result<Clinvar, Die> {
    let mut out = Clinvar {
        findings: Vec::new(),
        total: 0,
        aside: Vec::new(),
    };
    each_record(path, |r| {
        if r.len() < 8 {
            return Ok(());
        }
        out.total += 1;
        let owned: Vec<String> = r.iter().map(|s| (*s).to_owned()).collect();
        let h =
            Hit::parse(&owned, Floor::default()).map_err(|e| Die::one(format!("{path}: {e}")))?;
        let gt = r.get(9).map_or("", |s| s.split(':').next().unwrap_or(""));
        let pos = i64::try_from(h.pos)
            .map_err(|_| Die::one(format!("{path}: POS {} is out of range", h.pos)))?;
        let Some(z) = zygosity(gt, &h.chrom, sex, par.contains(&h.chrom, pos)) else {
            return Ok(());
        };
        let variant = format!("{}:{} {}>{}", h.chrom, h.pos, h.reference, h.alt);
        for g in h.gene.split(',') {
            let Some(gene) = list.get(g) else { continue };
            if let Some(sev) = gene.alleles.get(&h.key()) {
                let kind = sev.kind.as_deref().unwrap_or("");
                out.aside.push(format!(
                    "{variant} {} in {g}: {}, a {kind} allele",
                    z.word(),
                    sev.name
                ));
                if kind == "hypomorph" {
                    out.findings.push(Finding {
                        gene: g.to_owned(),
                        text: format!("{variant} {}, {}, a hypomorph allele", z.word(), sev.name),
                        copies: z.copies(),
                        hypomorph: Some(sev.name.clone()),
                    });
                }
                continue;
            }
            let ident = match r[2] {
                "." | "" => String::new(),
                id => format!(", ClinVar {id}"),
            };
            // A split half call counts as het (one copy) but reads as the
            // other reports word it.
            let words = if z == Zygosity::Het && half_call(gt) {
                format!("one copy, other copy not called ref ({gt})")
            } else {
                z.word().to_owned()
            };
            let text = format!(
                "{variant} {words}, {}, {}{ident}",
                h.significance,
                stars_text(h.stars)
            );
            out.findings.push(Finding::new(g, text, z.copies()));
        }
        Ok(())
    })?;
    Ok(out)
}

/// GENCODE exons as a BED (0-based start, end), merged per chromosome:
/// starts and ends, sorted.
pub struct Exons(HashMap<String, (Vec<i64>, Vec<i64>)>);

impl Exons {
    pub fn read(path: &str) -> Result<Exons, Die> {
        let mut by: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        each_line(path, Utf8::Replace, |line| {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 3 || line.starts_with('#') || line.starts_with("track") {
                return Ok(());
            }
            let s = int(c[1].trim(), path)?;
            let e = int(c[2].trim(), path)?;
            by.entry(c[0].to_owned()).or_default().push((s, e));
            Ok(())
        })?;
        let mut merged = HashMap::new();
        for (chrom, mut iv) in by {
            iv.sort_unstable();
            let mut m: Vec<(i64, i64)> = Vec::new();
            for (s, e) in iv {
                match m.last_mut() {
                    Some(last) if s <= last.1 => last.1 = last.1.max(e),
                    _ => m.push((s, e)),
                }
            }
            merged.insert(chrom, m.into_iter().unzip());
        }
        Ok(Exons(merged))
    }

    /// How many merged exons overlap the 0-based half-open `[s0, e0)`.
    pub fn hit(&self, chrom: &str, s0: i64, e0: i64) -> usize {
        let Some((starts, ends)) = self.0.get(chrom) else {
            return 0;
        };
        let i = ends.partition_point(|&e| e <= s0);
        starts[i..].iter().take_while(|&&s| s < e0).count()
    }
}

fn plain_chrom(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

/// AnnotSV's benign loss set: per chromosome (without `chr`), the records
/// of [`COMMON_AF`] or more, sorted.
pub struct Benign {
    by: HashMap<String, Vec<(i64, i64)>>,
    /// How many records are common.
    pub count: usize,
}

impl Benign {
    /// Reads the BED (chrom, start, end, source, coordinates, AF). A record
    /// without a numeric frequency is not taken as common.
    pub fn read(path: &str) -> Result<Benign, Die> {
        let mut by: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        let mut count = 0;
        each_line(path, Utf8::Replace, |line| {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 6
                || line.starts_with('#')
                || line.starts_with("track")
                || !is_digits(c[1])
                || !is_digits(c[2])
            {
                return Ok(());
            }
            let Ok(af) = pycompat::num::parse_float(c[5]) else {
                return Ok(());
            };
            if af >= COMMON_AF {
                by.entry(plain_chrom(c[0]).to_owned())
                    .or_default()
                    .push((int(c[1], path)?, int(c[2], path)?));
                count += 1;
            }
            Ok(())
        })?;
        for v in by.values_mut() {
            v.sort_unstable();
        }
        Ok(Benign { by, count })
    }

    /// Whether the 0-based half-open `[s0, e0)` and a common deletion each
    /// cover at least half of the other.
    pub fn matches(&self, chrom: &str, s0: i64, e0: i64) -> bool {
        let Some(rows) = self.by.get(plain_chrom(chrom)) else {
            return false;
        };
        let size = e0 - s0;
        // A match is at most twice as long, so it starts no earlier than that.
        let from = rows.partition_point(|r| r.0 < s0 - size * 100 / COMMON_OVERLAP_PERCENT - 1);
        rows[from..]
            .iter()
            .take_while(|r| r.0 < e0)
            .any(|&(bs, be)| {
                let ov = e0.min(be) - s0.max(bs);
                ov > 0
                    && ov * 100 >= COMMON_OVERLAP_PERCENT * size
                    && ov * 100 >= COMMON_OVERLAP_PERCENT * (be - bs)
            })
    }
}

/// The Delly deletions set aside over an exon of a listed gene.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct DellySkips {
    /// Shorter than [`SV_MIN_SIZE`].
    pub short: usize,
    /// Longer than [`DEPTH_MIN_SIZE`], and read depth does not show the loss.
    pub no_loss: usize,
    /// Matching a common deletion of the benign set.
    pub common: usize,
}

/// `re.search(r"(?:^|;)KEY=(VALUE)", info)`: the first `KEY=` segment
/// whose value has a match for `value` (the matched part).
fn info_first<'a>(info: &'a str, key: &str, value: fn(&str) -> &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| {
        let v = kv.strip_prefix(key)?.strip_prefix('=')?;
        let m = &v[..value(v).len()];
        (!m.is_empty()).then_some(m)
    })
}

/// `[^;]+`: a segment's value is all of it.
fn whole(v: &str) -> &str {
    v
}

/// `\d+`: the leading digits.
fn digits(v: &str) -> &str {
    &v[..v.bytes().take_while(u8::is_ascii_digit).count()]
}

fn sv_sex(sex: Option<&str>) -> Sex {
    match sex {
        Some("male") => Sex::Male,
        Some("female") => Sex::Female,
        _ => Sex::Unknown,
    }
}

/// PASS Delly deletions of [`SV_MIN_SIZE`] or more that remove at least
/// one exon of a listed gene, and the counts of those set aside.
pub fn delly(
    path: &str,
    index: &GeneIndex,
    exons: &Exons,
    sex: Option<&str>,
    par: &Par,
    benign: Option<&Benign>,
) -> Result<(Vec<Finding>, DellySkips), Die> {
    let mut out = Vec::new();
    let mut skip = DellySkips::default();
    let over_exon = |chrom: &str, pos: i64, end: i64| {
        index
            .overlapping(chrom, pos + 1, end)
            .iter()
            .any(|&(gs, ge, _)| exons.hit(chrom, pos.max(gs - 1), end.min(ge)) > 0)
    };
    each_record(path, |r| {
        if r.len() < 10 || r[6] != "PASS" {
            return Ok(());
        }
        let svtype = info_first(r[7], "SVTYPE", whole);
        let end = info_first(r[7], "END", digits);
        let (Some("DEL"), Some(end)) = (svtype, end) else {
            return Ok(());
        };
        let fmt: Vec<&str> = r[8].split(':').collect();
        let vals: Vec<&str> = r[9].split(':').collect();
        if let Some(i) = fmt.iter().position(|k| *k == "FT")
            && i < vals.len()
            && !matches!(vals[i], "PASS" | ".")
        {
            return Ok(());
        }
        let gt = if fmt[0] == "GT" { vals[0] } else { "" };
        let chrom = r[0];
        let pos = int(r[1], path)?;
        let end = int(end, path)?;
        let in_par = par.contains(chrom, pos);
        let Some(z) = zygosity(gt, chrom, sex, in_par) else {
            return Ok(());
        };
        let size = end - pos;
        if size < SV_MIN_SIZE {
            skip.short += usize::from(over_exon(chrom, pos, end));
            return Ok(());
        }
        // Read depth must show the loss of a long deletion: RDCN below the
        // copies a normal region has (og_steps::sv::depth_contradicts).
        if depth_contradicts(
            SvKind::Del,
            chrom,
            size,
            Rdcn::from_format(&fmt, &vals),
            sv_sex(sex),
            in_par,
        ) {
            skip.no_loss += usize::from(over_exon(chrom, pos, end));
            return Ok(());
        }
        // The deleted bases POS+1 to END are the 0-based [POS, END).
        if benign.is_some_and(|b| b.matches(chrom, pos, end)) {
            skip.common += usize::from(over_exon(chrom, pos, end));
            return Ok(());
        }
        for (gs, ge, g) in index.overlapping(chrom, pos + 1, end) {
            let n = exons.hit(chrom, pos.max(gs - 1), end.min(ge));
            if n == 0 {
                continue;
            }
            let span = if size >= 1000 {
                format!("{} kb", pycompat::fmt::fixed(size as f64 / 1000.0, 1))
            } else {
                format!("{size} bp")
            };
            let text = format!(
                "deletion {chrom}:{}-{end} ({span}) {}, {n} exon{} of the gene (Delly)",
                pos + 1,
                z.word(),
                if n == 1 { "" } else { "s" }
            );
            out.push(Finding::new(g, text, z.copies()));
        }
        Ok(())
    })?;
    Ok((out, skip))
}

/// FMR1 and FXN alleles in their expanded range (ExpansionHunter).
pub fn repeats(path: &str, list: &CarrierList) -> Result<Vec<Finding>, Die> {
    let mut out = Vec::new();
    each_record(path, |r| {
        if r.len() < 10 {
            return Ok(());
        }
        let locus = info_get(r[7], "VARID")
            .filter(|v| !v.is_empty())
            .or_else(|| info_get(r[7], "REPID"))
            .unwrap_or("");
        let Some(&(_, least, unit)) = REPEATS.iter().find(|(l, _, _)| *l == locus) else {
            return Ok(());
        };
        if !list.contains(locus) {
            return Ok(());
        }
        let fmt: Vec<&str> = r[8].split(':').collect();
        let vals: Vec<&str> = r[9].split(':').collect();
        let Some(cn) = fmt
            .iter()
            .position(|k| *k == "REPCN")
            .and_then(|i| vals.get(i))
        else {
            return Ok(());
        };
        let big: Vec<i64> = cn
            .split('/')
            .filter(|x| is_digits(x))
            .filter_map(|x| x.parse().ok())
            .filter(|&x| x >= least)
            .collect();
        let Some(&max) = big.iter().max() else {
            return Ok(());
        };
        let what = match locus {
            "FMR1" if max >= FMR1_FULL => "full mutation (more than 200)",
            "FMR1" => "premutation (55 to 200)",
            _ => "expanded (66 or more)",
        };
        let two = big.len() > 1;
        let text = format!(
            "{locus} {unit} repeat {cn}, {} allele{} {what} (ExpansionHunter estimate)",
            if two { "two" } else { "one" },
            if two { "s" } else { "" }
        );
        out.push(Finding::new(
            locus,
            text,
            u8::try_from(big.len()).unwrap_or(u8::MAX),
        ));
        Ok(())
    })?;
    Ok(out)
}

/// SMN1 copy number from step 35's table (the one rule, `og_cn::smn`):
/// the findings, and the text of what was read.
pub fn smn(path: &str, list: &CarrierList) -> Result<(Vec<Finding>, String), Die> {
    let p = og_cn::smn::read(Some(path), og_cn::smn::SMN1_EXON7.0, &[])
        .map_err(|e| Die::one(format!("{path}: {e}")))?;
    // Step 35 writes a row it could not assess with the filter "not
    // assessed: <reason>"; the reason is what the reader needs.
    if let Some(reason) = p
        .text
        .strip_suffix("), not resolved")
        .and_then(|t| t.split_once("filter not assessed: "))
        .map(|(_, r)| r)
    {
        return Ok((
            Vec::new(),
            format!("SMN1 not assessed: {reason}; copy number not resolved"),
        ));
    }
    let cn = match p.pscn {
        Some(cn) if list.contains("SMN1") => cn,
        _ => return Ok((Vec::new(), format!("SMN1 copy number {}", p.text))),
    };
    if cn >= 2 {
        return Ok((Vec::new(), format!("SMN1 {cn} copies")));
    }
    let word = if cn == 1 { "copy" } else { "copies" };
    let text = format!("SMN1 {cn} {word} (Parascopy: {})", p.text);
    Ok((
        vec![Finding::new("SMN1", text, if cn == 1 { 1 } else { 2 })],
        format!("SMN1 {cn} {word}"),
    ))
}

/// Step 35's CYP21A2 table (metric, value), a later metric over an earlier.
fn cyp21a2_table(path: &str) -> Result<HashMap<String, String>, Die> {
    let mut t = HashMap::new();
    each_line(path, Utf8::Replace, |line| {
        if let Some((k, v)) = line.split_once('\t')
            && k != "metric"
        {
            t.insert(k.to_owned(), v.to_owned());
        }
        Ok(())
    })?;
    Ok(t)
}

/// Python's `repr()` of an optional string, for a message.
fn py_repr(v: Option<&str>) -> String {
    match v {
        None => "None".to_owned(),
        Some(s) if s.contains('\'') && !s.contains('"') => format!("\"{s}\""),
        Some(s) => format!("'{}'", s.replace('\\', "\\\\").replace('\'', "\\'")),
    }
}

/// What step 35's CYP21A2 table says.
pub struct Cyp21a2 {
    pub findings: Vec<Finding>,
    /// Why it was not assessed, for the not-assessed list.
    pub gap: Option<String>,
}

/// CYP21A2 from depth and PSV sites (step 35's `<id>_cyp21a2_psv.tsv`). A
/// deletion or a possible conversion is a finding on one copy, on two when
/// a site has no CYP21A2 copy left.
pub fn cyp21a2(path: &str, list: &CarrierList) -> Result<Cyp21a2, Die> {
    let t = cyp21a2_table(path)?;
    let get = |k: &str| t.get(k).map_or(".", String::as_str);
    let status = t.get("status").map(String::as_str);
    if !matches!(status, Some("called" | "not assessed")) {
        return Err(Die::one(format!(
            "{path}: not step 35's CYP21A2 table (og cn call --locus CYP21A2)"
        )));
    }
    if status != Some("called") {
        let text = format!("CYP21A2 not assessed ({})", get("reason"));
        // The table's message says, for a BAM without multi-mapping reads,
        // that the aligner and not the genome stopped the call; an older
        // table has only the reason code.
        let gap = match get("message").strip_prefix("not assessed: ") {
            Some(_) => format!("CYP21A2 {}", get("message")),
            None => text.clone(),
        };
        return Ok(Cyp21a2 {
            findings: Vec::new(),
            gap: Some(gap),
        });
    }
    let what = format!(
        "CYP21A2 {} copies, RCCX module total {} of {}",
        get("cyp21a2_copies"),
        get("module_total"),
        get("expected_total")
    );
    let kind = t.get("finding").map(String::as_str);
    if !matches!(
        kind,
        Some("none" | "deletion" | "conversion" | "pseudogene_deletion")
    ) {
        return Err(Die::one(format!(
            "{path}: a called CYP21A2 table with the finding {} (none, deletion, conversion or pseudogene_deletion)",
            py_repr(kind)
        )));
    }
    let mut findings = Vec::new();
    if matches!(kind, Some("deletion" | "conversion")) && list.contains("CYP21A2") {
        let text = if kind == Some("deletion") {
            format!(
                "CYP21A2 deletion from depth and PSV sites ({what}): the 30 kb chimeric deletion pattern"
            )
        } else {
            format!(
                "possible CYP21A2 conversion at {} ({what}, PSV sites)",
                get("conversion_sites").replace("PSV_", "")
            )
        };
        let two = t.get("cyp21a2_copies").map(String::as_str) == Some("0");
        findings.push(Finding::new(
            "CYP21A2",
            format!("{text}; unconfirmed, short reads"),
            if two { 2 } else { 1 },
        ));
    }
    Ok(Cyp21a2 {
        findings,
        gap: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zygosity_counts_a_missing_allele_in_the_ploidy() {
        let z = |gt, chrom, sex| zygosity(gt, chrom, sex, false);
        // A split half call (./1, the GT 1/2 record step 06 splits) is one copy.
        assert_eq!(z("./1", "chr7", Some("female")), Some(Zygosity::Het));
        assert_eq!(z("1|.", "chr7", Some("male")), Some(Zygosity::Het));
        assert_eq!(z("1/1", "chr7", Some("female")), Some(Zygosity::Hom));
        assert_eq!(z("1", "chrX", Some("female")), Some(Zygosity::Hemizygous));
        assert_eq!(z("./1", "chrX", Some("male")), Some(Zygosity::Hemizygous));
        assert_eq!(z("./.", "chr7", Some("female")), None);
        assert_eq!(z("./0", "chr7", Some("female")), None);
        assert_eq!(z("", "chr7", None), None);
        assert!(half_call("./1") && half_call("2|.") && !half_call("0/1") && !half_call("./."));
    }

    #[test]
    fn exons_merge_and_count() {
        let mut m = HashMap::new();
        m.insert("chr1".to_owned(), (vec![10, 50], vec![20, 60]));
        let e = Exons(m);
        assert_eq!(e.hit("chr1", 0, 100), 2);
        assert_eq!(e.hit("chr1", 20, 50), 0);
        assert_eq!(e.hit("chr1", 19, 51), 2);
        assert_eq!(e.hit("chr2", 0, 100), 0);
    }

    #[test]
    fn info_first_takes_the_first_key() {
        assert_eq!(
            info_first("IMPRECISE;SVTYPE=DEL;END=5;END=9", "END", digits),
            Some("5")
        );
        assert_eq!(info_first("END=x;END=7z", "END", digits), Some("7"));
        assert_eq!(info_first("CIEND=1;XEND=3", "END", digits), None);
        assert_eq!(
            info_first("SVTYPE=;SVTYPE=DEL", "SVTYPE", whole),
            Some("DEL")
        );
        assert_eq!(info_get("A=1;A=2", "A"), Some("2"));
    }
}
