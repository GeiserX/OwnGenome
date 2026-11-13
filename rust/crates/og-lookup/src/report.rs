//! The answer for one region: a site's genotype, or a gene's variants and
//! how much of it has a call; then the annotation and the reports.

use crate::Error;
use crate::query::{contig_name, is_rsid, norm_chrom};
use crate::resolve::Row;
use crate::summary::{self, Summary};
use crate::vcf::{Csq, Entry, Rec, Region};
use pycompat::{fmt, text};
use std::collections::BTreeSet;

/// CSQ fields worth printing, in order, with their label.
const CSQ_SHOW: [(&str, &str); 12] = [
    ("Consequence", "consequence"),
    ("IMPACT", "impact"),
    ("HGVSc", "HGVSc"),
    ("HGVSp", "HGVSp"),
    ("Existing_variation", "known as"),
    ("ClinVar_CLNSIG", "ClinVar (pipeline's ClinVar file)"),
    ("ClinVar_CLNREVSTAT", "ClinVar review status"),
    ("ClinVar_CLNDN", "ClinVar condition"),
    ("CLIN_SIG", "clinical significance (VEP cache)"),
    ("gnomADg_AF", "gnomAD genomes AF"),
    ("gnomADe_AF", "gnomAD exomes AF"),
    ("MAX_AF", "highest population AF"),
];

/// What one region's answer reads.
pub struct Ask<'a> {
    /// The query as given.
    pub query: &'a str,
    pub row: &'a Row,
    pub gvcf: &'a Region,
    pub ann: Option<&'a Region>,
    pub summary: Option<&'a Summary>,
}

/// Appends the answer to `out`. On an error, `out` holds what was written
/// before it, as the Python had printed it.
pub fn write(a: &Ask, out: &mut String) -> Result<(), Error> {
    if a.row.kind == "site" {
        site(a, out)
    } else {
        gene(a, out)
    }
}

fn push(out: &mut String, line: impl AsRef<str>) {
    out.push_str(line.as_ref());
    out.push('\n');
}

/// The CSQ entry to print: of `gene` when given; then the most severe
/// IMPACT, and among those VEP's PICK, MANE Select, the canonical
/// transcript, else the first.
fn pick<'a>(entries: Vec<Entry<'a>>, gene: Option<&str>) -> Option<Entry<'a>> {
    let mut entries = entries;
    if let Some(g) = gene {
        let up = text::upper(g);
        let mine: Vec<Entry> = entries
            .iter()
            .filter(|e| text::upper(e.get("SYMBOL")) == up)
            .cloned()
            .collect();
        if !mine.is_empty() {
            entries = mine;
        }
    }
    let rank = |e: &Entry| match e.get("IMPACT") {
        "HIGH" => 0,
        "MODERATE" => 1,
        "LOW" => 2,
        "MODIFIER" => 3,
        _ => 4,
    };
    let best = entries.iter().map(rank).min()?;
    let entries: Vec<Entry> = entries.into_iter().filter(|e| rank(e) == best).collect();
    let tests: [&dyn Fn(&Entry) -> bool; 3] = [
        &|e| e.get("PICK") == "1",
        &|e| !e.get("MANE_SELECT").is_empty(),
        &|e| e.get("CANONICAL") == "YES",
    ];
    for t in tests {
        if let Some(e) = entries.iter().find(|e| t(e)) {
            return Some(e.clone());
        }
    }
    entries.into_iter().next()
}

/// heterozygous, homozygous, hemizygous, half call or no call.
fn zygosity<T: PartialEq>(bases: &[Option<T>]) -> &'static str {
    let called: Vec<&T> = bases.iter().flatten().collect();
    if bases.is_empty() || called.is_empty() {
        "no call"
    } else if called.len() < bases.len() {
        "half call"
    } else if bases.len() == 1 {
        "hemizygous"
    } else if called.iter().all(|b| *b == called[0]) {
        "homozygous"
    } else {
        "heterozygous"
    }
}

fn fmt_gt(bases: &[Option<&str>], phased: bool) -> String {
    if bases.is_empty() {
        return ".".into();
    }
    let parts: Vec<&str> = bases.iter().map(|b| b.unwrap_or(".")).collect();
    parts.join(if phased { "|" } else { "/" })
}

fn depth_text(r: &Rec) -> String {
    let mut parts = Vec::new();
    for (k, label) in [
        ("DP", "depth"),
        ("MIN_DP", "lowest depth in the block"),
        ("GQ", "genotype quality"),
    ] {
        if let Some(v) = r.fmt(k).filter(|v| !v.is_empty() && *v != ".") {
            parts.push(format!("{label} {v}"));
        }
    }
    if !r.qual.is_empty() && r.qual != "." && !r.is_block() {
        parts.push(format!("QUAL {}", r.qual));
    }
    if parts.is_empty() {
        "no depth or quality recorded".into()
    } else {
        parts.join(", ")
    }
}

fn filtered(r: &Rec) -> bool {
    !matches!(r.filter.as_str(), "PASS" | "." | "")
}

/// The genotype lines of one position of the gVCF.
fn site_genotype(gvcf: &Region, chrom: &str, pos: u64, ref_: &str) -> Result<Vec<String>, Error> {
    let cover: Vec<&Rec> = gvcf
        .records
        .iter()
        .filter(|r| r.chrom == chrom && r.pos <= pos && pos <= r.end)
        .collect();
    if cover.is_empty() {
        return Ok(vec![
            "genotype: no call (no gVCF record covers this position: no reads there, or the caller made no call)"
                .into(),
        ]);
    }
    let r = cover
        .iter()
        .find(|r| r.pos == pos && !r.is_block())
        .or_else(|| cover.iter().find(|r| r.is_block()))
        .unwrap_or(&cover[0]);
    let (bases, phased) = r.bases()?;
    let mut lines = Vec::new();
    let at = format!("{}:{}-{}", contig_name(&r.chrom), r.pos, r.end);
    if r.is_block() {
        let (idx, _) = r.gt()?;
        if idx.iter().all(Option::is_none) {
            lines.push(format!(
                "genotype: no call (reference block {at} has no genotype)"
            ));
        } else {
            // The table's base, else the block's own when it starts here;
            // an empty one reads as none.
            let base = if !ref_.is_empty() && ref_ != "." {
                ref_
            } else if r.pos == pos {
                r.ref_.as_str()
            } else {
                ""
            };
            let base = if base.is_empty() { "ref" } else { base };
            let is_ref = |i: &Option<i128>| *i == Some(0);
            let gt = idx
                .iter()
                .map(|i| if is_ref(i) { base } else { "." })
                .collect::<Vec<_>>()
                .join("/");
            let what = if idx.iter().all(is_ref) {
                "homozygous reference"
            } else {
                let marks: Vec<Option<u8>> = idx.iter().map(|i| is_ref(i).then_some(0)).collect();
                zygosity(&marks)
            };
            lines.push(format!(
                "genotype: {gt} ({what}, inside the reference block {at})"
            ));
        }
    } else if r.pos != pos {
        lines.push(format!(
            "genotype: this position lies inside the variant called at {}:{} {}>{}, genotype {} ({})",
            contig_name(&r.chrom),
            r.pos,
            r.ref_,
            r.real_alts().join(","),
            fmt_gt(&bases, phased),
            zygosity(&bases)
        ));
    } else {
        let alts = r.real_alts().join(",");
        lines.push(format!(
            "genotype: {} ({}; REF {}, ALT {})",
            fmt_gt(&bases, phased),
            zygosity(&bases),
            r.ref_,
            if alts.is_empty() { "." } else { &alts }
        ));
    }
    lines.push(format!("evidence: {}", depth_text(r)));
    if filtered(r) {
        lines.push(format!(
            "filter: {} (the caller did not pass this call; read it as uncertain)",
            r.filter
        ));
    }
    Ok(lines)
}

/// VEP writes a comma inside a value as `&`; ClinVar writes spaces as `_`.
/// Then `re.sub(r"\s*,\s*", ", ", v).strip()`.
fn clinvar_text(v: &str) -> String {
    let v = v.replace('&', ",").replace('_', " ");
    let parts: Vec<&str> = v.split(',').collect();
    let last = parts.len() - 1;
    let parts: Vec<&str> = parts
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let p = if i > 0 {
                p.trim_start_matches(text::is_space)
            } else {
                p
            };
            if i < last {
                p.trim_end_matches(text::is_space)
            } else {
                p
            }
        })
        .collect();
    text::strip(&parts.join(", "), None).to_owned()
}

fn annotation_lines(r: &Rec, fields: &Csq, gene: Option<&str>, indent: &str) -> Vec<String> {
    let e = pick(r.csq(fields), gene);
    let mut out = Vec::new();
    let get = |k: &str| e.as_ref().map_or("", |e| e.get(k));
    let sym = get("SYMBOL");
    if !sym.is_empty() {
        out.push(format!("{indent}gene: {sym}"));
    }
    for (k, label) in CSQ_SHOW {
        let v = get(k);
        let v = if k.starts_with("ClinVar_") || k == "CLIN_SIG" {
            clinvar_text(v)
        } else {
            v.replace('&', ", ")
        };
        if !v.is_empty() {
            out.push(format!("{indent}{label}: {v}"));
        }
    }
    for k in ["CLNSIG", "CLNREVSTAT", "CLNDN"] {
        if let Some(v) = r.info_str(k).filter(|v| !v.is_empty())
            && get(&format!("ClinVar_{k}")).is_empty()
        {
            out.push(format!("{indent}ClinVar {k}: {}", v.replace('_', " ")));
        }
    }
    if !fields.is_empty() && e.is_none() {
        out.push(format!("{indent}(no VEP annotation on this record)"));
    }
    out
}

fn site(a: &Ask, out: &mut String) -> Result<(), Error> {
    let chrom = norm_chrom(&a.row.chrom);
    let pos = a.row.start;
    let from = if a.row.source.is_empty() {
        String::new()
    } else {
        format!("; from {}", a.row.source)
    };
    push(
        out,
        format!(
            "== {}: {}:{pos} (GRCh38{from})",
            a.query,
            contig_name(&chrom)
        ),
    );
    for line in site_genotype(a.gvcf, &chrom, pos, &a.row.ref_)? {
        push(out, format!("  {line}"));
    }
    let mut genes = BTreeSet::new();
    let mut rsid = is_rsid(a.query).then(|| a.query.to_ascii_lowercase());
    match a.ann {
        None => push(
            out,
            "  annotation: no annotated VCF (step 13), so no consequence, ClinVar or gnomAD",
        ),
        Some(ann) => {
            let recs: Vec<&Rec> = ann
                .records
                .iter()
                .filter(|r| r.chrom == chrom && r.pos <= pos && pos <= r.end)
                .collect();
            if recs.is_empty() {
                push(
                    out,
                    "  annotation: no variant of this sample here in the annotated VCF (it holds only the \
                     sample's variant calls); ClinVar and gnomAD are not read for a reference or no-call site",
                );
            }
            for r in recs {
                let (bases, phased) = r.bases()?;
                push(
                    out,
                    format!(
                        "  variant: {}:{} {}>{} {} ({})",
                        contig_name(&r.chrom),
                        r.pos,
                        r.ref_,
                        r.alts.join(","),
                        fmt_gt(&bases, phased),
                        zygosity(&bases)
                    ),
                );
                for line in annotation_lines(r, &ann.csq, None, "    ") {
                    push(out, line);
                }
                genes.extend(r.symbols(&ann.csq));
                if rsid.is_none() {
                    rsid = r
                        .ids(&ann.csq)
                        .iter()
                        .map(|i| text::lower(i))
                        .find(|i| i.starts_with("rs"));
                }
            }
        }
    }
    push(out, "  in the reports:");
    for line in summary::lines(a.summary, &genes, rsid.as_deref())? {
        push(out, line);
    }
    Ok(())
}

/// The bases of `[start, end]` that some called gVCF record covers: the
/// union of the records' spans, so overlapping records (a deletion and the
/// reference block after it) count each base once.
fn covered(gvcf: &Region, chrom: &str, start: u64, end: u64) -> Result<u64, Error> {
    let mut spans = Vec::new();
    for r in &gvcf.records {
        if r.chrom != chrom || r.end < start || r.pos > end {
            continue;
        }
        let (idx, _) = r.gt()?;
        if idx.iter().any(Option::is_some) {
            let (b, e) = (r.pos.max(start), r.end.min(end));
            if b <= e {
                spans.push((b, e));
            }
        }
    }
    spans.sort_unstable();
    let mut total = 0;
    let mut reach = 0; // the last base counted so far
    for (b, e) in spans {
        let b = b.max(reach + 1);
        if b <= e {
            total += e - b + 1;
            reach = e;
        }
    }
    Ok(total)
}

fn gene(a: &Ask, out: &mut String) -> Result<(), Error> {
    let chrom = norm_chrom(&a.row.chrom);
    let (start, end) = (a.row.start, a.row.end);
    let gene = a.query;
    let source = if a.row.source.is_empty() {
        "gene span"
    } else {
        &a.row.source
    };
    push(
        out,
        format!(
            "== {gene}: {}:{start}-{end} (GRCh38; {source})",
            contig_name(&chrom)
        ),
    );
    let span = (end + 1)
        .checked_sub(start)
        .filter(|&s| s > 0)
        .ok_or_else(|| Error::other(format!("the gene span {start}-{end} is empty")))?;
    let cov = covered(a.gvcf, &chrom, start, end)?;
    push(
        out,
        format!(
            "  gVCF: {} of {} bases have a genotype call ({}%); a variant can only be missing from the list \
             below where there is no call",
            fmt::int_grouped(i128::from(cov)),
            fmt::int_grouped(i128::from(span)),
            fmt::fixed(100.0 * cov as f64 / span as f64, 1)
        ),
    );
    let mut listed = 0;
    let inside = |r: &&Rec| r.chrom == chrom && r.end >= start && r.pos <= end;
    if let Some(ann) = a.ann {
        push(out, "  variants called in the gene (annotated VCF):");
        let up = text::upper(gene);
        for r in ann.records.iter().filter(inside) {
            if !ann.csq.is_empty() && !r.symbols(&ann.csq).iter().any(|s| text::upper(s) == up) {
                continue;
            }
            let (bases, phased) = r.bases()?;
            let ids = r.ids(&ann.csq).join(",");
            let flag = if filtered(r) {
                format!(", filter {}", r.filter)
            } else {
                String::new()
            };
            push(
                out,
                format!(
                    "  - {}:{} {}>{} {} ({}; {}; {}{flag})",
                    contig_name(&r.chrom),
                    r.pos,
                    r.ref_,
                    r.alts.join(","),
                    fmt_gt(&bases, phased),
                    zygosity(&bases),
                    if ids.is_empty() { "no rsID" } else { &ids },
                    depth_text(r)
                ),
            );
            for line in annotation_lines(r, &ann.csq, Some(gene), "      ") {
                push(out, line);
            }
            listed += 1;
        }
    } else {
        push(
            out,
            "  variants called in the gene (gVCF; no annotated VCF from step 13, so no consequence, ClinVar or gnomAD):",
        );
        for r in a.gvcf.records.iter().filter(inside) {
            if r.is_block() || !r.has_alt_call()? {
                continue;
            }
            let (bases, phased) = r.bases()?;
            let flag = if filtered(r) {
                format!(", filter {}", r.filter)
            } else {
                String::new()
            };
            push(
                out,
                format!(
                    "  - {}:{} {}>{} {} ({}; {}{flag})",
                    contig_name(&r.chrom),
                    r.pos,
                    r.ref_,
                    r.real_alts().join(","),
                    fmt_gt(&bases, phased),
                    zygosity(&bases),
                    depth_text(r)
                ),
            );
            listed += 1;
        }
    }
    if listed == 0 {
        push(
            out,
            "  - none: the sample matches the reference wherever the gene has a call",
        );
    }
    push(out, "  in the reports:");
    let genes = BTreeSet::from([gene.to_owned()]);
    for line in summary::lines(a.summary, &genes, None)? {
        push(out, line);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clinvar_text_is_python_re_sub() {
        assert_eq!(
            clinvar_text("criteria_provided&_multiple_submitters"),
            "criteria provided, multiple submitters"
        );
        assert_eq!(clinvar_text("a , b,c "), "a, b, c");
        assert_eq!(clinvar_text("a, ,b"), "a, , b");
        assert_eq!(clinvar_text(""), "");
    }

    #[test]
    fn coverage_is_a_union() {
        let rec = |l: &str| Rec::parse(&l.split('\t').collect::<Vec<_>>()).unwrap();
        let g = Region {
            csq: Csq::default(),
            records: vec![
                rec("chr7\t100\t.\tACGTA\tA\t50\tPASS\t.\tGT\t0/1"),
                rec("chr7\t101\t.\tC\t<*>\t0\t.\tEND=103\tGT\t0/0"),
                rec("chr7\t106\t.\tC\t<*>\t0\t.\tEND=106\tGT\t./."),
            ],
        };
        // 100-104 by the deletion, 101-103 again by the block: 5 bases, not 8.
        assert_eq!(covered(&g, "7", 100, 109).unwrap(), 5);
    }

    #[test]
    fn zygosity_names() {
        assert_eq!(zygosity(&[Some("A"), Some("T")]), "heterozygous");
        assert_eq!(zygosity(&[Some("A"), None]), "half call");
        assert_eq!(zygosity(&[Some("A")]), "hemizygous");
        assert_eq!(zygosity::<&str>(&[None, None]), "no call");
        assert_eq!(zygosity::<&str>(&[]), "no call");
    }
}
