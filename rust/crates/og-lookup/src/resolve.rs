//! From a query to the regions to read: rsIDs through the consumer site
//! table, the pipeline's tables and the annotated VCF; genes through
//! GENCODE's gene lines.
//!
//! Each source is read once for every query of a command line: the site
//! table and the annotated VCF are each one pass, whatever the number of
//! rsIDs.

use crate::query::{Query, contig_name, norm_chrom};
use crate::vcf::{Csq, Rec, int};
use crate::{Error, Inputs};
use memchr::memmem;
use og_io::{Lines, Utf8};
use pycompat::text;
use std::collections::HashMap;
use std::io::BufRead;

/// One region to read for a query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    /// `site` or `gene`.
    pub kind: &'static str,
    /// The normalised chromosome.
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// The GRCh38 base when a table gives it, else `.`.
    pub ref_: String,
    /// Where the region came from.
    pub source: String,
}

impl Row {
    /// The line `og lookup resolve` prints.
    pub fn line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            self.kind,
            contig_name(&self.chrom),
            self.start,
            self.end,
            self.ref_,
            self.source
        )
    }
}

/// A position an rsID source names: chromosome, position, base.
type Hit = (String, u64, String);

/// Hits per rsID, or why the source could not be read.
type Hits = Result<HashMap<String, Vec<Hit>>, Error>;

/// The annotated VCF's hits per rsID, each with its own error: a record
/// that names one rsID and cannot be parsed stops only that rsID's query,
/// as the Python's one scan per query did.
type PerRsid = Result<HashMap<String, Result<Vec<Hit>, Error>>, Error>;

const SITES_SOURCE: &str = "the consumer site table (dbSNP positions)";
const ANNOTATED_SOURCE: &str = "the annotated VCF (a variant called in this sample)";

/// Resolves every query of one command line, in order: each one's rows, or
/// the error that stops it (exit 3: found nowhere; 2: an unreadable query
/// or table).
pub fn resolve_all(queries: &[String], inp: &Inputs) -> Vec<Result<Vec<Row>, Error>> {
    let parsed: Vec<Result<Query, Error>> =
        queries.iter().map(|q| crate::query::parse(q)).collect();
    let rsids: Vec<&str> = parsed
        .iter()
        .filter_map(|q| match q {
            Ok(Query::Rsid(r)) => Some(r.as_str()),
            _ => None,
        })
        .collect();
    let genes: Vec<&str> = parsed
        .iter()
        .filter_map(|q| match q {
            Ok(Query::Gene(g)) => Some(g.as_str()),
            _ => None,
        })
        .collect();

    let sites = match &inp.sites {
        Some(p) if !rsids.is_empty() => Some(sites_rsid(p, &rsids)),
        _ => None,
    };
    let tables: Vec<(String, Hits)> = if rsids.is_empty() {
        Vec::new()
    } else {
        inp.tables
            .iter()
            .map(|t| {
                let name = t.rsplit('/').next().unwrap_or(t).to_owned();
                (name, table_rsid(t, &rsids))
            })
            .collect()
    };
    // The annotated VCF is read only for the rsIDs no table placed.
    let found_in_tables = |r: &str| -> bool {
        let hit = |h: &Hits| {
            h.as_ref()
                .is_ok_and(|m| m.get(r).is_some_and(|v| !v.is_empty()))
        };
        sites.as_ref().is_some_and(hit) || tables.iter().any(|(_, h)| hit(h))
    };
    let unplaced: Vec<&str> = rsids
        .iter()
        .copied()
        .filter(|r| !found_in_tables(r))
        .collect();
    let annotated = match &inp.annotated {
        Some(p) if !unplaced.is_empty() => Some(annotated_rsid(p, &unplaced)),
        _ => None,
    };
    let spans = match &inp.genes {
        Some(p) if !genes.is_empty() => Some(gene_spans(p, &genes)),
        _ => None,
    };

    parsed
        .into_iter()
        .map(|q| match q? {
            Query::Site(c, p) => Ok(vec![Row {
                kind: "site",
                chrom: c,
                start: p,
                end: p,
                ref_: ".".into(),
                source: "position as given".into(),
            }]),
            Query::Rsid(r) => rsid_rows(&r, inp, &sites, &tables, &annotated),
            Query::Gene(g) => gene_rows(&g, inp, &spans),
        })
        .collect()
}

fn rsid_rows(
    rsid: &str,
    inp: &Inputs,
    sites: &Option<Hits>,
    tables: &[(String, Hits)],
    annotated: &Option<PerRsid>,
) -> Result<Vec<Row>, Error> {
    let take = |h: &Hits| -> Result<Vec<Hit>, Error> {
        match h {
            Ok(m) => Ok(m.get(rsid).cloned().unwrap_or_default()),
            Err(e) => Err(e.clone()),
        }
    };
    let mut hits: Vec<(Hit, String)> = Vec::new();
    if let Some(s) = sites {
        hits.extend(take(s)?.into_iter().map(|h| (h, SITES_SOURCE.to_owned())));
    }
    for (name, t) in tables {
        hits.extend(take(t)?.into_iter().map(|h| (h, name.clone())));
    }
    if hits.is_empty()
        && let Some(a) = annotated
    {
        let found = match a {
            Ok(m) => m.get(rsid).cloned().unwrap_or(Ok(Vec::new()))?,
            Err(e) => return Err(e.clone()),
        };
        hits.extend(found.into_iter().map(|h| (h, ANNOTATED_SOURCE.to_owned())));
    }
    let mut seen = Vec::new();
    let mut rows = Vec::new();
    for ((c, p, base), source) in hits {
        if seen.contains(&(c.clone(), p)) {
            continue;
        }
        seen.push((c.clone(), p));
        rows.push(Row {
            kind: "site",
            chrom: c,
            start: p,
            end: p,
            ref_: if base.is_empty() { ".".into() } else { base },
            source,
        });
    }
    if rows.is_empty() {
        let mut where_ = vec![if inp.sites.is_some() {
            "the consumer site table"
        } else {
            "no consumer site table (setup.sh --consumer-sites builds one)"
        }];
        if !inp.tables.is_empty() {
            where_.push("the pipeline's tables");
        }
        where_.push(if inp.annotated.is_some() {
            "the annotated VCF"
        } else {
            "no annotated VCF (step 13)"
        });
        return Err(Error::not_found(format!(
            "{rsid} was not found in {}. An rsID the sample does not carry is only known through the site \
             table; look it up on dbSNP (https://www.ncbi.nlm.nih.gov/snp/) and ask for its GRCh38 position \
             instead, as chr:pos.",
            where_.join(", ")
        )));
    }
    Ok(rows)
}

fn gene_rows(
    gene: &str,
    inp: &Inputs,
    spans: &Option<Result<Spans, Error>>,
) -> Result<Vec<Row>, Error> {
    let (Some(path), Some(spans)) = (&inp.genes, spans) else {
        return Err(Error::not_found(
            "no GENCODE gene coordinates (reference/gencode.v*.basic.genes.gtf): run scripts/setup.sh, \
             which installs them, or ask for a position (chr:pos).",
        ));
    };
    let spans = spans.as_ref().map_err(Clone::clone)?;
    let rows = spans.for_gene(gene);
    if rows.is_empty() {
        return Err(Error::not_found(format!(
            "{gene} is not a gene name in {path}. Use the HGNC symbol (for example CFTR, not its alias), an \
             rsID or a position (chr:pos)."
        )));
    }
    Ok(rows)
}

/// The lines of a small text file as Python's `open()` yields them, each
/// without its `\n`.
fn text_lines(path: &str) -> Result<impl Iterator<Item = Result<String, Error>>, Error> {
    let f = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let p = path.to_owned();
    Ok(Lines::new(f, Utf8::Replace).map(move |l| l.map_err(|e| Error::io(&e, &p))))
}

/// The consumer site table's rows for `rsids`. It opens with `## key:
/// value` lines and a `#rsid ...` column line.
fn sites_rsid(path: &str, rsids: &[&str]) -> Hits {
    let mut out: HashMap<String, Vec<Hit>> = HashMap::new();
    let mut head: Option<HashMap<String, usize>> = None;
    for line in text_lines(path)? {
        let line = line?;
        if line.starts_with("##") || text::strip(&line, None).is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let Some(h) = &head else {
            let h: HashMap<String, usize> = c
                .iter()
                .enumerate()
                .map(|(i, k)| (k.trim_start_matches('#').to_owned(), i))
                .collect();
            if !["rsid", "chrom", "pos38", "ref38"]
                .iter()
                .all(|k| h.contains_key(*k))
            {
                return Err(Error::value(format!(
                    "{path}: not a consumer site table (columns {})",
                    c.join(", ")
                )));
            }
            head = Some(h);
            continue;
        };
        let col = |k: &str| {
            c.get(h[k])
                .copied()
                .ok_or_else(|| Error::other(format!("{path}: a row has no {k} column")))
        };
        let id = text::lower(col("rsid")?);
        if rsids.contains(&id.as_str()) {
            let pos = to_pos(int(col("pos38")?)?)?;
            out.entry(id).or_default().push((
                norm_chrom(col("chrom")?),
                pos,
                col("ref38")?.to_owned(),
            ));
        }
    }
    Ok(out)
}

/// A pipeline table's rows for `rsids`: `assets/traits.tsv` (chrom pos
/// rsid ref) or `assets/bloodtype_markers.tsv` (a BED: its end is the
/// 1-based position). Lines starting with `#` are comments; the first
/// other line is the header.
fn table_rsid(path: &str, rsids: &[&str]) -> Hits {
    let mut out: HashMap<String, Vec<Hit>> = HashMap::new();
    let mut head: Option<Vec<String>> = None;
    for line in text_lines(path)? {
        let line = line?;
        if line.starts_with('#') || text::strip(&line, None).is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let Some(h) = &head else {
            head = Some(c.iter().map(|s| (*s).to_owned()).collect());
            continue;
        };
        // dict(zip(head, c)): the last column of a repeated name wins.
        let get = |k: &str| -> Option<&str> {
            let n = h.len().min(c.len());
            (0..n).rev().find(|&i| h[i] == k).map(|i| c[i])
        };
        let id = text::lower(get("rsid").unwrap_or(""));
        if !rsids.contains(&id.as_str()) {
            continue;
        }
        let pos = match get("pos").or_else(|| get("end")) {
            Some(p) => to_pos(int(p)?)?,
            None => continue,
        };
        let r = get("ref").filter(|r| !r.is_empty()).unwrap_or(".");
        let base = if r == "." {
            ".".to_owned()
        } else {
            r.chars().take(1).collect()
        };
        let chrom = get("chrom").ok_or_else(|| Error::other(format!("{path}: no chrom column")))?;
        out.entry(id)
            .or_default()
            .push((norm_chrom(chrom), pos, base));
    }
    Ok(out)
}

fn to_pos(v: i128) -> Result<u64, Error> {
    u64::try_from(v).map_err(|_| Error::value(format!("a negative position: {v}")))
}

/// The annotated VCF's records whose ID or Existing_variation is one of
/// `rsids`, in one pass over the whole file (BGZF inflated on several
/// cores). A line is parsed only when it holds one of the rsIDs' digits
/// right after `rs` in any case; the record's IDs then decide.
fn annotated_rsid(path: &str, rsids: &[&str]) -> PerRsid {
    let mut out: HashMap<String, Result<Vec<Hit>, Error>> = HashMap::new();
    let finders: Vec<memmem::Finder> = rsids
        .iter()
        .map(|r| memmem::Finder::new(&r.as_bytes()[2..]))
        .collect();
    let mut file =
        og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(|e| Error::io(&e, path))?;
    let mut head: Vec<String> = Vec::new();
    let mut fields = Csq::default();
    let mut raw = Vec::new();
    loop {
        raw.clear();
        if file
            .read_until(b'\n', &mut raw)
            .map_err(|e| Error::io(&e, path))?
            == 0
        {
            break;
        }
        let body = raw.strip_suffix(b"\n").unwrap_or(&raw);
        let body = body.strip_suffix(b"\r").unwrap_or(body);
        // Most lines name none of the rsIDs and are never decoded. A lone
        // \r ends a line in Python's text mode, so such a line is split.
        let has_cr = memchr::memchr(b'\r', body).is_some();
        if body.first() != Some(&b'#') && !has_cr && !names_one(&finders, body) {
            continue;
        }
        let text = String::from_utf8_lossy(body);
        for line in text.split('\r') {
            if line.starts_with('#') {
                head.push(line.to_owned());
                if line.starts_with("#CHROM") {
                    fields = Csq::from_header(head.iter().map(String::as_str));
                }
                continue;
            }
            // The rsIDs this line names: a line the Python's scan for an
            // rsID matched, and so parsed, for that rsID's query only.
            let named: Vec<&str> = rsids
                .iter()
                .zip(&finders)
                .filter(|(_, f)| names(f, line.as_bytes()))
                .map(|(r, _)| *r)
                .collect();
            if named.is_empty() {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            let r = match Rec::parse(&c) {
                Ok(r) => r,
                Err(e) => {
                    for rsid in named {
                        let slot = out.entry(rsid.to_owned()).or_insert_with(|| Ok(Vec::new()));
                        if slot.is_ok() {
                            *slot = Err(e.clone());
                        }
                    }
                    continue;
                }
            };
            let ids: Vec<String> = r.ids(&fields).iter().map(|i| text::lower(i)).collect();
            for rsid in named {
                if ids.iter().any(|i| i == rsid)
                    && let Ok(v) = out.entry(rsid.to_owned()).or_insert_with(|| Ok(Vec::new()))
                {
                    let base: String = r.ref_.chars().take(1).collect();
                    v.push((r.chrom.clone(), r.pos, base));
                }
            }
        }
    }
    Ok(out)
}

/// Whether `line` names one of the rsIDs as a word: `rs` in any case and
/// the digits, with no letter or digit before and no digit after (the
/// boundary rule the Python's lookaround regex had).
fn names_one(finders: &[memmem::Finder], line: &[u8]) -> bool {
    finders.iter().any(|f| names(f, line))
}

/// Whether `line` names the one rsID `f` finds, by the rule of [`names_one`].
fn names(f: &memmem::Finder, line: &[u8]) -> bool {
    let n = f.needle().len();
    f.find_iter(line).any(|at| {
        at >= 2
            && line[at - 2..at].eq_ignore_ascii_case(b"rs")
            && (at == 2 || !line[at - 3].is_ascii_alphanumeric())
            && !line.get(at + n).is_some_and(u8::is_ascii_digit)
    })
}

/// GENCODE's gene lines, by name.
struct Spans {
    source: String,
    /// (name, chrom, start, end) in file order.
    genes: Vec<(String, String, u64, u64)>,
}

impl Spans {
    /// The spans of `gene`: its exact name, else a case-insensitive match.
    fn for_gene(&self, gene: &str) -> Vec<Row> {
        let row = |(n, c, s, e): &(String, String, u64, u64)| Row {
            kind: "gene",
            chrom: c.clone(),
            start: *s,
            end: *e,
            ref_: ".".into(),
            source: format!("{} ({n})", self.source),
        };
        let exact: Vec<Row> = self.genes.iter().filter(|g| g.0 == gene).map(row).collect();
        if !exact.is_empty() {
            return exact;
        }
        let up = text::upper(gene);
        self.genes
            .iter()
            .filter(|g| text::upper(&g.0) == up)
            .map(row)
            .collect()
    }
}

/// The gene lines of `path` named like one of `genes` (case-insensitive).
fn gene_spans(path: &str, genes: &[&str]) -> Result<Spans, Error> {
    let source = match path.find("gencode.v") {
        Some(at) => {
            let d: String = path[at + 9..]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect();
            if d.is_empty() {
                "GENCODE".to_owned()
            } else {
                format!("GENCODE {d}")
            }
        }
        None => "GENCODE".to_owned(),
    };
    let wanted: Vec<String> = genes.iter().map(|g| text::upper(g)).collect();
    let mut out = Vec::new();
    for line in text_lines(path)? {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 9 || c[2] != "gene" {
            continue;
        }
        let Some(name) = gene_name(c[8]) else {
            continue;
        };
        if !wanted.contains(&text::upper(name)) {
            continue;
        }
        out.push((
            name.to_owned(),
            norm_chrom(c[0]),
            to_pos(int(c[3])?)?,
            to_pos(int(c[4])?)?,
        ));
    }
    Ok(Spans { source, genes: out })
}

/// The value of `gene_name "..."` (`gene_name "([^"]+)"`).
fn gene_name(attrs: &str) -> Option<&str> {
    let mut rest = attrs;
    while let Some(at) = rest.find("gene_name \"") {
        let tail = &rest[at + 11..];
        match tail.find('"') {
            Some(q) if q > 0 => return Some(&tail[..q]),
            _ => rest = &rest[at + 1..],
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_rsid_is_found_as_a_word_only() {
        let f = [memmem::Finder::new(b"1001")];
        assert!(names_one(&f, b"chr7\t1\trs1001\tA"));
        assert!(names_one(&f, b"x|RS1001&COSV1|y"));
        assert!(!names_one(&f, b"chr7\t1\trs10011\tA"));
        assert!(!names_one(&f, b"chr7\t1\tcrs1001\tA"));
        assert!(!names_one(&f, b"chr7\t1\t1001\tA"));
    }

    #[test]
    fn a_bad_record_stops_only_the_rsid_it_names() {
        let path =
            std::env::temp_dir().join(format!("og-lookup-badpos-{}.vcf", std::process::id()));
        std::fs::write(
            &path,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr7\tx\trs1\tA\tG\t.\tPASS\t.\nchr7\t20\trs2\tC\tT\t.\tPASS\t.\n",
        )
        .unwrap();
        let inp = Inputs {
            gvcf: None,
            annotated: Some(path.to_string_lossy().into_owned()),
            summary: None,
            sites: None,
            tables: Vec::new(),
            genes: None,
        };
        let got = resolve_all(&["rs1".into(), "rs2".into()], &inp);
        std::fs::remove_file(&path).unwrap();
        assert!(got[0].is_err(), "rs1 names the bad line: {:?}", got[0]);
        let rows = got[1].as_ref().expect("rs2 is not on the bad line");
        assert_eq!((rows.len(), rows[0].start), (1, 20));
    }

    #[test]
    fn gene_name_takes_the_first_quoted_value() {
        assert_eq!(
            gene_name("gene_id \"E\"; gene_name \"CFTR\"; level 2;"),
            Some("CFTR")
        );
        assert_eq!(gene_name("gene_name \"\"; gene_name \"X\";"), Some("X"));
        assert_eq!(gene_name("gene_id \"E\";"), None);
    }
}
