//! `og dev build acmg-lof`: the `clingen_hi` and `lof` columns of
//! `assets/acmg_sf.tsv`, and TTN's `limit` cell, built from two pinned
//! downloads and one reviewed file.
//!
//! `scripts/ci/build-acmg-lof.sh` fetches the pinned downloads (archived
//! copies, so their bytes never change) and checks their sha256; this verb
//! only reads them. It rewrites those cells and leaves every other byte of
//! the table as it is, so the table is its own oracle: `--check` exits 1
//! when a rebuild would change it.
//!
//! `clingen_hi` is ClinGen's haploinsufficiency score as ClinGen writes it
//! (`3` sufficient evidence, `2` some, `1` little, `0` none, `30` a gene of
//! an autosomal recessive disease, `40` dosage sensitivity unlikely), or `.`
//! for a gene ClinGen has not curated. `lof` says whether a predicted
//! loss-of-function call (a VEP HIGH-impact record) counts as a secondary
//! finding in the gene (`og summary collect`, `collect/acmg.rs`). A score
//! decides it alone only where it proves it:
//!
//! - `3` (sufficient evidence that losing one copy causes disease) is `yes`;
//! - `40` (dosage sensitivity unlikely) is `no`: a PCSK9 stop-gain is no
//!   finding, as PCSK9 causes familial hypercholesterolemia by gain of
//!   function.
//!
//! Every other score (`0`, `1`, `2`, `30`, `.`) says too little, so the gene
//! needs a row of `assets/acmg_lof_overrides.tsv` (its `lof` and the
//! reason), reviewed like code; a gene with neither stops the build. A row
//! there also overrides a score of 3 or 40 that belongs to another disease
//! (APOB, RET). So a gene a later ACMG list adds can never lose its
//! loss-of-function findings without someone writing down why.
//!
//! TTN's `limit` cell is `truncating_in` the TTN exons with high PSI: ACMG
//! SF v3.2, Table 1 footnote f, counts only "frameshift and nonsense
//! variants, and variants known to impact the splicing of TTN exons with
//! high PSI", citing Roberts et al. 2015, whose high is PSI > 0.9. The PSI
//! of each exon is cardiodb's titin exon table (the PSI of each exon in
//! the left ventricles of 84 DCM patients, Roberts et al. 2015): each run of
//! neighbouring exons above 90% becomes one region, widened by 2 bases at
//! each end for the canonical splice sites.

use super::{Args, BuildError, input_error};
use indexmap::IndexMap;
use std::io::Write;
use std::path::Path;

const USAGE: &str = "usage: og dev build acmg-lof --clingen FILE --overrides FILE --ttn-exons FILE --table FILE [--check]\n";

/// The gene whose `limit` cell the TTN exon table gives.
const TTN: &str = "TTN";

/// The name the `limit` cell gives TTN's regions.
const TTN_NAME: &str = "the TTN exons with PSI over 90%";

/// An exon counts when more than this percentage of TTN transcripts in the
/// heart use it (Roberts et al. 2015: PSI > 0.9).
const PSI_OVER: u32 = 90;

/// The exons of cardiodb's table: the 364 exons of the TTN gene in LRG
/// numbering (LRG_391), one more than ENST00000589042 has (exon 48 is
/// novex-3's last exon).
const TTN_EXONS: u32 = 364;

/// Bases on each side of a run of exons that belong to its canonical
/// splice sites.
const SPLICE: i128 = 2;

/// The `lof` cell of a gene: the hand row of `overrides` when it has one,
/// else what its ClinGen score proves, else `None` (it needs a hand row).
pub fn lof<'a>(gene: &str, hi: &str, overrides: &'a Overrides) -> Option<&'a str> {
    if let Some((l, _)) = overrides.get(gene) {
        return Some(l.as_str());
    }
    match hi {
        "3" => Some("yes"),
        "40" => Some("no"),
        _ => None,
    }
}

/// `{gene: (lof, reason)}` of `assets/acmg_lof_overrides.tsv`.
pub type Overrides = IndexMap<String, (String, String)>;

/// Reads `assets/acmg_lof_overrides.tsv`: `#` lines are comments, then a
/// header naming `gene`, `lof` and `reason`, then one row per gene.
pub fn overrides(path: &Path) -> Result<Overrides, BuildError> {
    let text = std::fs::read_to_string(path).map_err(|e| input_error(path, &e))?;
    let name = path.display();
    let mut cols: Option<[usize; 3]> = None;
    let mut out = Overrides::new();
    for line in text.lines() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let Some([g, l, r]) = cols else {
            let at = |c: &str| {
                f.iter()
                    .position(|x| *x == c)
                    .ok_or_else(|| BuildError::Input(format!("{name}: no {c} column")))
            };
            cols = Some([at("gene")?, at("lof")?, at("reason")?]);
            continue;
        };
        let cell = |i: usize| f.get(i).map_or("", |x| x.trim());
        let (gene, l, reason) = (cell(g), cell(l), cell(r));
        if gene.is_empty() {
            return Err(BuildError::Input(format!(
                "{name}: a row has no gene: {line:?}"
            )));
        }
        if !matches!(l, "yes" | "no" | "mark") {
            return Err(BuildError::Data(format!(
                "{name}: {gene}: lof must be yes, no or mark, not '{l}'"
            )));
        }
        if matches!(reason, "" | ".") {
            return Err(BuildError::Data(format!("{name}: {gene} has no reason")));
        }
        if out
            .insert(gene.to_owned(), (l.to_owned(), reason.to_owned()))
            .is_some()
        {
            return Err(BuildError::Data(format!("{name}: {gene} is on two rows")));
        }
    }
    if cols.is_none() {
        return Err(BuildError::Input(format!("{name}: no header row")));
    }
    Ok(out)
}

/// `{symbol: haploinsufficiency score}` of ClinGen's gene curation list:
/// `#` lines are comments, except the one that names the columns.
pub fn clingen_scores(path: &Path) -> Result<IndexMap<String, String>, BuildError> {
    let text = std::fs::read_to_string(path).map_err(|e| input_error(path, &e))?;
    let name = path.display();
    let mut cols: Option<(usize, usize)> = None;
    let mut out = IndexMap::new();
    for line in text.lines() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if let Some(h) = line.strip_prefix("#Gene Symbol") {
            let names: Vec<&str> = std::iter::once("Gene Symbol")
                .chain(h.split('\t').skip(1))
                .collect();
            let at = |c: &str| names.iter().position(|n| *n == c);
            cols = Some((
                0,
                at("Haploinsufficiency Score").ok_or_else(|| {
                    BuildError::Input(format!("{name}: no Haploinsufficiency Score column"))
                })?,
            ));
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let (g, s) = cols.ok_or_else(|| {
            BuildError::Input(format!(
                "{name}: a gene row comes before the #Gene Symbol line"
            ))
        })?;
        let f: Vec<&str> = line.split('\t').collect();
        let (Some(gene), Some(score)) = (f.get(g), f.get(s)) else {
            return Err(BuildError::Input(format!(
                "{name}: a row has too few columns: {line:?}"
            )));
        };
        let score = score.trim();
        if gene.is_empty() {
            return Err(BuildError::Input(format!(
                "{name}: a row has no gene symbol"
            )));
        }
        let score = if score.is_empty() { "." } else { score };
        if out.insert((*gene).to_owned(), score.to_owned()).is_some() {
            return Err(BuildError::Data(format!("{name}: {gene} is on two rows")));
        }
    }
    if cols.is_none() {
        return Err(BuildError::Input(format!("{name}: no #Gene Symbol line")));
    }
    Ok(out)
}

/// One row of cardiodb's titin exon table.
#[derive(Debug, Clone, PartialEq)]
pub struct Exon {
    /// The LRG exon number.
    pub number: u32,
    /// Its GRCh38 ends on chr2, in either order (TTN is on the minus
    /// strand, so the table writes the higher one first).
    pub ends: (i128, i128),
    /// Its PSI in percent; `None` when the table gives none (exon 48,
    /// "n/a").
    pub psi: Option<u32>,
}

/// The text of an HTML cell: tags dropped, `&nbsp;` and `&amp;` read,
/// blanks collapsed.
fn cell_text(html: &str) -> String {
    let mut out = String::new();
    let mut in_tag = false;
    for c in html.chars() {
        match c {
            '<' => in_tag = true,
            '>' if in_tag => in_tag = false,
            _ if !in_tag => out.push(c),
            _ => {}
        }
    }
    out.replace("&nbsp;", " ")
        .replace("&amp;", "&")
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// The inner HTML of each `<tag ...>...</tag>` of `row`, in order.
fn cells<'a>(row: &'a str, tag: &str) -> Vec<&'a str> {
    let (open, close) = (format!("<{tag}"), format!("</{tag}>"));
    let mut out = Vec::new();
    let mut rest = row;
    while let Some(i) = rest.find(&open) {
        let after = &rest[i + open.len()..];
        // `<th` must not match `<thead`.
        if !after.starts_with(['>', ' ', '\t', '\n', '\r']) {
            rest = after;
            continue;
        }
        let Some(gt) = after.find('>') else { break };
        let body = &after[gt + 1..];
        let Some(end) = body.find(&close) else { break };
        out.push(&body[..end]);
        rest = &body[end + close.len()..];
    }
    out
}

/// The columns of cardiodb's exon table this builder reads, as the first
/// header row names them: the exon number, the GRCh38 ends and the PSI.
/// The seven columns of "Exon Rank in Transcript" sit between the meta
/// coordinates and the PSI.
const HEADER: [&str; 10] = [
    "Exon",
    "Hg19 start",
    "Hg19 end",
    "Hg38",
    "Hg38 end",
    "Meta start",
    "Meta end",
    "Exon Rank in Transcript",
    "PSI%",
    "Band",
];

/// The exons of cardiodb's titin exon table (the table `ttn_overview` of
/// its page `titin_transcripts.php`), numbered 1 to 364 in order.
pub fn ttn_exons(html: &str) -> Result<Vec<Exon>, BuildError> {
    let bad = |m: String| BuildError::Input(format!("the TTN exon table: {m}"));
    let start = html
        .find("<table id=\"ttn_overview\"")
        .ok_or_else(|| bad("no table ttn_overview".into()))?;
    let table = &html[start..];
    let table = &table[..table
        .find("</table>")
        .ok_or_else(|| bad("no </table>".into()))?];
    let mut rows = table.split("<tr").skip(1);
    let head: Vec<String> = rows
        .next()
        .map(|r| cells(r, "th").into_iter().map(cell_text).collect())
        .unwrap_or_default();
    let ok = head.len() >= HEADER.len()
        && HEADER
            .iter()
            .zip(&head)
            .all(|(want, got)| got.starts_with(want));
    if !ok {
        return Err(bad(format!(
            "its header is not Exon, Hg19 start, Hg19 end, Hg38 start, Hg38 end, Meta start, \
             Meta end, Exon Rank in Transcript, PSI%, Band: {head:?}"
        )));
    }
    let mut out: Vec<Exon> = Vec::new();
    for row in rows {
        let c: Vec<String> = cells(row, "td").into_iter().map(cell_text).collect();
        // Exon number, 2 hg19, 2 hg38, 2 meta, 7 ranks, PSI, band.
        if c.len() < 16 {
            continue;
        }
        let Ok(number) = c[0].parse::<u32>() else {
            continue;
        };
        let pos = |i: usize| {
            c[i].parse::<i128>()
                .ok()
                .filter(|p| *p > 0)
                .ok_or_else(|| bad(format!("exon {number}: '{}' is no position", c[i])))
        };
        let psi = if c[14].chars().all(|x| x.is_ascii_digit()) && !c[14].is_empty() {
            let p = c[14]
                .parse::<u32>()
                .map_err(|_| bad(format!("exon {number}: PSI '{}'", c[14])))?;
            if p > 100 {
                return Err(bad(format!("exon {number}: PSI {p} is over 100")));
            }
            Some(p)
        } else {
            None
        };
        let want = out.len() as u32 + 1;
        if number != want {
            return Err(bad(format!("exon {number} comes where exon {want} should")));
        }
        out.push(Exon {
            number,
            ends: (pos(3)?, pos(4)?),
            psi,
        });
    }
    if out.len() as u32 != TTN_EXONS {
        return Err(bad(format!(
            "{} exons, not the {TTN_EXONS} of TTN",
            out.len()
        )));
    }
    Ok(out)
}

/// The GRCh38 chr2 regions of `exons` with PSI over 90%: one per run of
/// neighbouring exons, from its lowest to its highest base, widened by 2
/// bases at each end; sorted by position. An exon with no PSI breaks a run.
pub fn high_psi_regions(exons: &[Exon]) -> Vec<(i128, i128)> {
    let mut out: Vec<(i128, i128)> = Vec::new();
    let mut run: Option<(i128, i128)> = None;
    for e in exons {
        let (lo, hi) = (e.ends.0.min(e.ends.1), e.ends.0.max(e.ends.1));
        if e.psi.is_some_and(|p| p > PSI_OVER) {
            run = Some(run.map_or((lo, hi), |(a, b)| (a.min(lo), b.max(hi))));
        } else if let Some(r) = run.take() {
            out.push(r);
        }
    }
    out.extend(run);
    let mut out: Vec<(i128, i128)> = out
        .into_iter()
        .map(|(a, b)| (a - SPLICE, b + SPLICE))
        .collect();
    out.sort_unstable();
    out
}

/// TTN's `limit` cell from its exon table.
pub fn ttn_limit(exons: &[Exon]) -> Result<String, BuildError> {
    let regions = high_psi_regions(exons);
    if regions.is_empty() {
        return Err(BuildError::Data(
            "the TTN exon table has no exon with PSI over 90%".into(),
        ));
    }
    let list: Vec<String> = regions
        .iter()
        .map(|(a, b)| format!("chr2:{a}-{b}"))
        .collect();
    Ok(format!("truncating_in {} {TTN_NAME}", list.join(",")))
}

/// The table with its `clingen_hi`, `lof` and TTN `limit` cells rebuilt,
/// and the genes whose cells changed. A gene whose score settles nothing
/// and that has no hand row, or a hand row whose gene has no table row,
/// stops it.
pub fn rebuild(
    table: &str,
    scores: &IndexMap<String, String>,
    overrides: &Overrides,
    ttn: &str,
) -> Result<(String, Vec<String>), BuildError> {
    let mut out = String::with_capacity(table.len());
    let mut header: Option<[usize; 4]> = None;
    let mut changed = Vec::new();
    let mut seen = Vec::new();
    let mut missing = Vec::new();
    for line in table.split_inclusive('\n') {
        if line.starts_with('#') {
            out.push_str(line);
            continue;
        }
        let body = line.trim_end_matches(['\n', '\r']);
        let end = &line[body.len()..];
        let mut cells: Vec<String> = body.split('\t').map(str::to_owned).collect();
        let Some([g, hi, lf, lim]) = header else {
            let at = |c: &str| {
                cells
                    .iter()
                    .position(|x| x == c)
                    .ok_or_else(|| BuildError::Input(format!("acmg_sf.tsv has no {c} column")))
            };
            header = Some([at("gene")?, at("clingen_hi")?, at("lof")?, at("limit")?]);
            out.push_str(line);
            continue;
        };
        if body.is_empty() {
            out.push_str(line);
            continue;
        }
        if cells.len() <= [g, hi, lf, lim].into_iter().max().unwrap_or(0) {
            return Err(BuildError::Input(format!(
                "acmg_sf.tsv: a row has too few columns: {body:?}"
            )));
        }
        let gene = cells[g].clone();
        let score = scores.get(&gene).map_or(".", String::as_str).to_owned();
        let Some(l) = lof(&gene, &score, overrides).map(str::to_owned) else {
            missing.push(format!("{gene} (score {score})"));
            out.push_str(line);
            continue;
        };
        let limit = if gene == TTN {
            ttn.to_owned()
        } else {
            cells[lim].clone()
        };
        if cells[hi] != score || cells[lf] != l || cells[lim] != limit {
            changed.push(gene.clone());
        }
        cells[hi] = score;
        cells[lf] = l;
        cells[lim] = limit;
        seen.push(gene);
        out.push_str(&cells.join("\t"));
        out.push_str(end);
    }
    if header.is_none() {
        return Err(BuildError::Input("acmg_sf.tsv has no header row".into()));
    }
    if !missing.is_empty() {
        return Err(BuildError::Data(format!(
            "ClinGen's score settles only 3 (lof yes) and 40 (lof no); give each of these genes \
             a row with its lof and the reason in assets/acmg_lof_overrides.tsv: {}",
            missing.join(", ")
        )));
    }
    if !seen.iter().any(|s| s == TTN) {
        return Err(BuildError::Data(
            "acmg_sf.tsv has no TTN row for the TTN exon table".into(),
        ));
    }
    let stale: Vec<&str> = overrides
        .keys()
        .filter(|g| !seen.iter().any(|s| s == *g))
        .map(String::as_str)
        .collect();
    if !stale.is_empty() {
        return Err(BuildError::Data(format!(
            "assets/acmg_lof_overrides.tsv has a row for a gene with no row in acmg_sf.tsv: {}",
            stale.join(", ")
        )));
    }
    Ok((out, changed))
}

fn run(args: &[String], err: &mut dyn Write) -> Result<u8, Result<BuildError, String>> {
    let a = Args::parse(
        args,
        &["--clingen", "--overrides", "--ttn-exons", "--table"],
        &["--check"],
    )
    .map_err(Err)?;
    let need = |f: &str| a.need(f).map_err(Err);
    let (clingen, hand, exons, table) = (
        need("--clingen")?,
        need("--overrides")?,
        need("--ttn-exons")?,
        need("--table")?,
    );
    let scores = clingen_scores(Path::new(clingen)).map_err(Ok)?;
    let hand = overrides(Path::new(hand)).map_err(Ok)?;
    let exons_path = Path::new(exons);
    let html = std::fs::read(exons_path).map_err(|e| Ok(input_error(exons_path, &e)))?;
    let ttn = ttn_limit(&ttn_exons(&String::from_utf8_lossy(&html)).map_err(Ok)?).map_err(Ok)?;
    let path = Path::new(table);
    let old = std::fs::read_to_string(path).map_err(|e| Ok(input_error(path, &e)))?;
    let (new, changed) = rebuild(&old, &scores, &hand, &ttn).map_err(Ok)?;
    if a.has("--check") {
        if new == old {
            return Ok(0);
        }
        let _ = writeln!(
            err,
            "{table}: not what the pinned downloads and assets/acmg_lof_overrides.tsv give: {} \
             (run scripts/ci/build-acmg-lof.sh)",
            if changed.is_empty() {
                "the rows' bytes differ".to_owned()
            } else {
                format!(
                    "clingen_hi, lof or limit differs for {}",
                    changed.join(", ")
                )
            }
        );
        return Ok(1);
    }
    super::write_atomic(path, new.as_bytes()).map_err(Ok)?;
    let _ = writeln!(err, "{table}: {} gene row(s) changed", changed.len());
    Ok(0)
}

/// `og dev build acmg-lof ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if matches!(args.first().map(String::as_str), Some("-h" | "--help")) {
        let _ = out.write_all(USAGE.as_bytes());
        return 0;
    }
    match run(args, err) {
        Ok(code) => code,
        Err(Err(m)) => {
            let _ = writeln!(err, "{USAGE}og dev build acmg-lof: error: {m}");
            2
        }
        Err(Ok(e)) => {
            let _ = writeln!(err, "ERROR: {e}");
            e.code()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CLINGEN: &str = "#ClinGen Gene Curation Results\n#07 Sep,2026\n\
#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\n\
BRCA2\t675\t13q13.1\tchr13:1-2\t3\tSufficient evidence for dosage pathogenicity\n\
PCSK9\t255738\t1p32.3\tchr1:1-2\t40\tDosage sensitivity unlikely\n\
MUTYH\t4595\t1p34.1\tchr1:3-4\t30\tGene associated with autosomal recessive phenotype\n\
MYH7\t4625\t14q11.2\tchr14:1-2\t0\tNo evidence available\n\
APOB\t338\t2p24.1\tchr2:1-2\t3\tSufficient evidence\n\
DSG2\t1829\t18q12.1\tchr18:1-2\t2\tSome evidence\n\
TTN\t7273\t2q31.2\tchr2:3-4\t2\tSome evidence\n";

    const HAND: &str = "# hand rows\ngene\tlof\treason\n\
MUTYH\tyes\trecessive: both copies lost\n\
MYH7\tno\tdominant-negative missense\n\
APOB\tno\tits score is for another disease\n\
DSG2\tyes\ttruncating variants cause ARVC\n\
TTN\tyes\tfootnote f\n\
GAA\tyes\trecessive: both copies lost\n";

    const TABLE: &str = "# a comment\ngene\tinheritance\tnote\tclingen_hi\tlof\tlimit\n\
BRCA2\tAD\t.\t\t\t.\nPCSK9\tAD\t.\t\t\t.\nMUTYH\tAR\t.\t\t\t.\nMYH7\tAD\t.\t\t\t.\n\
APOB\tAD\t.\t\t\t.\nDSG2\tAD\t.\t\t\t.\nTTN\tAD\t.\t\t\t.\nGAA\tAR\t.\t\t\t.\n";

    const TTN_CELL: &str = "truncating_in chr2:8-22 the TTN exons with PSI over 90%";

    /// A fresh folder per call: the tests run in parallel in one process.
    fn folder(tag: &str) -> std::path::PathBuf {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let d = std::env::temp_dir().join(format!("og-acmg-lof-{tag}-{}-{n}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn read<T>(text: &str, f: impl Fn(&Path) -> Result<T, BuildError>) -> Result<T, BuildError> {
        let d = folder("in");
        let p = d.join("in.tsv");
        std::fs::write(&p, text).unwrap();
        let out = f(&p);
        std::fs::remove_dir_all(&d).unwrap();
        out
    }

    fn scores() -> IndexMap<String, String> {
        read(CLINGEN, clingen_scores).unwrap()
    }

    fn hand(text: &str) -> Overrides {
        read(text, overrides).unwrap()
    }

    fn cells_of(table: &str, gene: &str) -> (String, String, String) {
        let row = table
            .lines()
            .find(|l| l.starts_with(&format!("{gene}\t")))
            .unwrap();
        let f: Vec<&str> = row.split('\t').collect();
        (f[3].to_owned(), f[4].to_owned(), f[5].to_owned())
    }

    #[test]
    fn a_score_settles_3_and_40_and_a_hand_row_everything_else() {
        let (t, changed) = rebuild(TABLE, &scores(), &hand(HAND), TTN_CELL).unwrap();
        let want = [
            ("BRCA2", "3", "yes"),
            ("PCSK9", "40", "no"),
            ("MUTYH", "30", "yes"),
            ("MYH7", "0", "no"),
            ("APOB", "3", "no"),
            ("DSG2", "2", "yes"),
            ("TTN", "2", "yes"),
            ("GAA", ".", "yes"),
        ];
        for (g, hi, l) in want {
            let (h, f, _) = cells_of(&t, g);
            assert_eq!((h.as_str(), f.as_str()), (hi, l), "{g}");
        }
        assert_eq!(changed.len(), want.len());
        // TTN's limit comes from the exon table; no other limit moves.
        assert_eq!(cells_of(&t, "TTN").2, TTN_CELL);
        assert_eq!(cells_of(&t, "BRCA2").2, ".");
        // The comment, the header and the other cells keep their bytes.
        assert!(t.starts_with("# a comment\ngene\tinheritance\tnote\tclingen_hi\tlof\tlimit\n"));
        assert!(t.contains("\nBRCA2\tAD\t.\t3\tyes\t.\n"));
    }

    /// A gene whose score is not 3 or 40, or that ClinGen has not curated,
    /// stops the build until it has a hand row: an AD gene a later ACMG
    /// list adds cannot silently lose its loss-of-function findings.
    #[test]
    fn a_gene_the_score_does_not_settle_needs_a_hand_row() {
        for (gene, score) in [("DES", "."), ("DSG2", "2"), ("MYH7", "0"), ("MUTYH", "30")] {
            let table = format!(
                "{TABLE}{}",
                if gene == "DES" {
                    "DES\tAD\t.\t\t\t.\n"
                } else {
                    ""
                }
            );
            let h: String = HAND
                .lines()
                .filter(|l| !l.starts_with(&format!("{gene}\t")))
                .map(|l| format!("{l}\n"))
                .collect();
            let e = rebuild(&table, &scores(), &hand(&h), TTN_CELL).unwrap_err();
            assert!(
                matches!(&e, BuildError::Data(m) if m.contains(&format!("{gene} (score {score})"))
                    && m.contains("acmg_lof_overrides.tsv")),
                "{gene}: {e}"
            );
        }
        // With its hand row the same gene builds, and the row decides.
        let table = format!("{TABLE}DES\tAD\t.\t\t\t.\n");
        let h = format!("{HAND}DES\tno\tdominant-negative missense\n");
        let (t, _) = rebuild(&table, &scores(), &hand(&h), TTN_CELL).unwrap();
        assert_eq!(cells_of(&t, "DES").1, "no");
        let h = format!("{HAND}DES\tmark\tless certain\n");
        let (t, _) = rebuild(&table, &scores(), &hand(&h), TTN_CELL).unwrap();
        assert_eq!(cells_of(&t, "DES").1, "mark");
        // A hand row overrides a score of 3 too (APOB), and only a hand row
        // turns a 3 into no.
        let h: String = HAND
            .lines()
            .filter(|l| !l.starts_with("APOB\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        let (t, _) = rebuild(TABLE, &scores(), &hand(&h), TTN_CELL).unwrap();
        assert_eq!(cells_of(&t, "APOB").1, "yes");
    }

    #[test]
    fn a_rebuild_is_byte_stable_and_check_sees_a_changed_cell() {
        let (s, h) = (scores(), hand(HAND));
        let (once, _) = rebuild(TABLE, &s, &h, TTN_CELL).unwrap();
        let (twice, changed) = rebuild(&once, &s, &h, TTN_CELL).unwrap();
        assert_eq!(once, twice);
        assert!(changed.is_empty());
        // CRLF endings stay CRLF.
        let crlf = once.replace('\n', "\r\n");
        assert_eq!(rebuild(&crlf, &s, &h, TTN_CELL).unwrap().0, crlf);
        let edited = once.replace("PCSK9\tAD\t.\t40\tno", "PCSK9\tAD\t.\t40\tyes");
        assert_eq!(rebuild(&edited, &s, &h, TTN_CELL).unwrap().1, ["PCSK9"]);
        let edited_ttn = once.replace("chr2:8-22", "chr2:8-23");
        assert_eq!(rebuild(&edited_ttn, &s, &h, TTN_CELL).unwrap().1, ["TTN"]);

        let d = folder("cli");
        std::fs::write(d.join("clingen.tsv"), CLINGEN).unwrap();
        std::fs::write(d.join("hand.tsv"), HAND).unwrap();
        std::fs::write(d.join("exons.html"), exon_page(&psi_rows())).unwrap();
        let ttn = ttn_limit(&ttn_exons(&exon_page(&psi_rows())).unwrap()).unwrap();
        let table = d.join("acmg_sf.tsv");
        let args = |check: bool| {
            let mut a: Vec<String> = ["--clingen", "--overrides", "--ttn-exons", "--table"]
                .iter()
                .zip([
                    d.join("clingen.tsv"),
                    d.join("hand.tsv"),
                    d.join("exons.html"),
                    table.clone(),
                ])
                .flat_map(|(f, p)| [(*f).to_owned(), p.display().to_string()])
                .collect();
            if check {
                a.push("--check".into());
            }
            a
        };
        let call = |check: bool| {
            let (mut o, mut e) = (Vec::new(), Vec::new());
            let code = cli(&args(check), &mut o, &mut e);
            (code, String::from_utf8(e).unwrap())
        };
        std::fs::write(&table, TABLE).unwrap();
        let (code, e) = call(true);
        assert_eq!(code, 1, "{e}");
        assert!(e.contains("PCSK9") && e.contains("TTN"), "{e}");
        assert_eq!(call(false).0, 0);
        let built = std::fs::read_to_string(&table).unwrap();
        assert_eq!(built, rebuild(TABLE, &s, &h, &ttn).unwrap().0);
        assert_eq!(call(true), (0, String::new()));
        std::fs::write(
            &table,
            built.replace("PCSK9\tAD\t.\t40\tno", "PCSK9\tAD\t.\t40\tyes"),
        )
        .unwrap();
        let (code, e) = call(true);
        assert_eq!(code, 1);
        assert!(e.contains("differs for PCSK9"), "{e}");
        // A gene with no hand row fails the build: exit 1, named.
        std::fs::write(
            d.join("hand.tsv"),
            HAND.replace("DSG2\tyes\ttruncating variants cause ARVC\n", ""),
        )
        .unwrap();
        let (code, e) = call(false);
        assert_eq!(code, 1, "{e}");
        assert!(e.contains("DSG2 (score 2)"), "{e}");
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn bad_inputs_are_refused() {
        let (s, h) = (scores(), hand(HAND));
        assert!(matches!(
            rebuild("gene\tclingen_hi\tlof\nA\t3\tyes\n", &s, &h, TTN_CELL),
            Err(BuildError::Input(m)) if m.contains("limit")
        ));
        // A hand row whose gene left the table is stale.
        let no_gaa: String = TABLE
            .lines()
            .filter(|l| !l.starts_with("GAA"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert!(
            matches!(rebuild(&no_gaa, &s, &h, TTN_CELL), Err(BuildError::Data(m)) if m.contains("GAA"))
        );
        let no_ttn = TABLE.replace("TTN\tAD\t.\t\t\t.\n", "");
        let h2: Overrides = h
            .iter()
            .filter(|(g, _)| *g != "TTN")
            .map(|(g, v)| (g.clone(), v.clone()))
            .collect();
        assert!(
            matches!(rebuild(&no_ttn, &s, &h2, TTN_CELL), Err(BuildError::Data(m)) if m.contains("TTN"))
        );
        // The hand rows: a known lof, a reason, one row per gene, a header.
        for bad in [
            "gene\tlof\treason\nA\tmaybe\twhy\n",
            "gene\tlof\treason\nA\tyes\t\n",
            "gene\tlof\treason\nA\tyes\t.\n",
            "gene\tlof\treason\nA\tyes\twhy\nA\tno\twhy\n",
        ] {
            assert!(
                matches!(read(bad, overrides), Err(BuildError::Data(_))),
                "{bad}"
            );
        }
        assert!(matches!(
            read("gene\tlof\n", overrides),
            Err(BuildError::Input(_))
        ));
        assert!(matches!(
            read("# only a comment\n", overrides),
            Err(BuildError::Input(_))
        ));
        let d = folder("bad");
        let p = d.join("c.tsv");
        std::fs::write(&p, format!("{CLINGEN}BRCA2\t1\tx\ty\t3\tz\n")).unwrap();
        assert!(matches!(clingen_scores(&p), Err(BuildError::Data(m)) if m.contains("BRCA2")));
        std::fs::write(&p, "BRCA2\t1\n").unwrap();
        assert!(matches!(clingen_scores(&p), Err(BuildError::Input(_))));
        std::fs::remove_dir_all(&d).unwrap();
        let (mut o, mut e) = (Vec::new(), Vec::new());
        assert_eq!(cli(&["--table".into(), "x".into()], &mut o, &mut e), 2);
        assert!(String::from_utf8(e).unwrap().contains("--clingen"));
    }

    /// 364 exons on the minus strand, 100 bases each with 100-base introns,
    /// exon 1 highest; PSI 100 except where `psi` says.
    fn psi_rows() -> Vec<(u32, &'static str)> {
        vec![
            (3, "90"),
            (4, "12"),
            (6, "n/a<sup>4</sup>"),
            (7, "96"),
            (364, "&nbsp;"),
        ]
    }

    fn exon_page(psi: &[(u32, &str)]) -> String {
        let mut rows = String::new();
        for n in 1..=TTN_EXONS {
            let hi = 1_000_000 - i128::from(n) * 200;
            let p = psi.iter().find(|(e, _)| *e == n).map_or("100", |(_, v)| v);
            let td = |x: &str| format!("<td><a href=\"titin_exon.php?id={n}\">{x}</a></td>");
            rows.push_str(&format!(
                "<tr>\n    {}{}{}{}{}{}{}{}{}{}{}{}\n    </tr>",
                td(&n.to_string()),
                td("1"),
                td("2"),
                td(&hi.to_string()),
                td(&(hi - 99).to_string()),
                td(""),
                td(""),
                (0..7).map(|_| td("1")).collect::<String>(),
                td(p),
                td("A-band"),
                td("<img src=\"x.png\">"),
                td("100"),
            ));
        }
        format!(
            "<html><table id=\"domain_icon\"><tr><th>x</th></tr></table>\
             <table id=\"ttn_overview\"><tr>\n<th rowspan=\"2\">Exon<span class=\"superscript\">\
             <a href=\"#footnotes\">1</a></span><br>Number</th>\n<th rowspan=\"2\">Hg19 start</th>\
             <th rowspan=\"2\">Hg19 end</th>\n<th rowspan=\"2\">Hg38<span>5</span> start</th>\
             <th rowspan=\"2\">Hg38 end</th>\n<th rowspan=\"2\">Meta start</th><th rowspan=\"2\">\
             Meta end</th>\n<th colspan=\"7\">Exon Rank in Transcript</th>\n<th rowspan=\"2\">\
             PSI%<span>2</span></th>\n<th rowspan=\"2\">Band</th>\n<th>Domains</th></tr>\n\
             <tr><th class=\"transcript_header\">Meta</th></tr>{rows}</table></html>"
        )
    }

    /// Exons 1-2 and 5 and 7-363 are over 90%; 3 (90, not over), 4 (12),
    /// 6 (n/a) and 364 (no PSI) break the runs. Each run spans its exons
    /// and 2 bases on each side.
    #[test]
    fn the_ttn_regions_are_the_runs_of_exons_over_90_percent() {
        let exons = ttn_exons(&exon_page(&psi_rows())).unwrap();
        assert_eq!(exons.len(), 364);
        assert_eq!(exons[0].ends, (999_800, 999_701));
        assert_eq!(exons[5].psi, None);
        assert_eq!(exons[2].psi, Some(90));
        let top = |n: i128| 1_000_000 - n * 200;
        let bottom = |n: i128| top(n) - 99;
        assert_eq!(
            high_psi_regions(&exons),
            [
                (bottom(363) - 2, top(7) + 2),
                (bottom(5) - 2, top(5) + 2),
                (bottom(2) - 2, top(1) + 2),
            ]
        );
        assert_eq!(
            ttn_limit(&exons).unwrap(),
            format!(
                "truncating_in chr2:{}-{},chr2:{}-{},chr2:{}-{} the TTN exons with PSI over 90%",
                bottom(363) - 2,
                top(7) + 2,
                bottom(5) - 2,
                top(5) + 2,
                bottom(2) - 2,
                top(1) + 2
            )
        );
        // No exon over 90%: nothing to count, so no limit.
        let none: Vec<(u32, &str)> = (1..=TTN_EXONS).map(|n| (n, "50")).collect();
        let low = ttn_exons(&exon_page(&none)).unwrap();
        assert!(matches!(ttn_limit(&low), Err(BuildError::Data(_))));
    }

    #[test]
    fn a_changed_exon_table_is_refused() {
        let page = exon_page(&psi_rows());
        // A column moved: the header no longer names PSI% where it was.
        let moved = page.replace("PSI%<span>2</span>", "Usage");
        assert!(matches!(ttn_exons(&moved), Err(BuildError::Input(m)) if m.contains("header")));
        // An exon missing, or out of order.
        let gap = page.replace("titin_exon.php?id=7\">7<", "titin_exon.php?id=7\">8<");
        assert!(matches!(ttn_exons(&gap), Err(BuildError::Input(m)) if m.contains("exon 8")));
        let short = page.replacen("<tr>\n    <td>", "<tr>\n    <td>x</td><td>", 0);
        assert!(ttn_exons(&short).is_ok());
        assert!(matches!(
            ttn_exons("<html></html>"),
            Err(BuildError::Input(_))
        ));
        let over = page.replace(">96<", ">960<");
        assert!(matches!(ttn_exons(&over), Err(BuildError::Input(m)) if m.contains("over 100")));
        let nopos = page.replacen("999800", "x", 1);
        assert!(
            matches!(ttn_exons(&nopos), Err(BuildError::Input(m)) if m.contains("no position"))
        );
        // Fewer than 364 exons.
        let cut = page.rsplit_once("<tr>").unwrap().0.to_owned() + "</table>";
        assert!(matches!(ttn_exons(&cut), Err(BuildError::Input(m)) if m.contains("363 exons")));
    }

    /// The shipped table agrees with the shipped hand rows: each gene's
    /// `lof` is what its `clingen_hi` and its hand row give, every gene the
    /// score does not settle has a hand row with a reason, and every hand
    /// row names a gene of the table. The full rebuild from the pinned
    /// downloads runs in the e2e case i411-sf-mechanism-1-acmg-lof.
    #[test]
    fn the_shipped_lof_column_follows_the_hand_rows() {
        let root = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets"));
        let hand = overrides(&root.join("acmg_lof_overrides.tsv")).unwrap();
        let table = std::fs::read_to_string(root.join("acmg_sf.tsv")).unwrap();
        let mut rows = table.lines().filter(|l| !l.starts_with('#'));
        let head: Vec<&str> = rows.next().unwrap().split('\t').collect();
        let at = |c: &str| head.iter().position(|x| *x == c).unwrap();
        let (g, hi, lf, lim) = (at("gene"), at("clingen_hi"), at("lof"), at("limit"));
        let mut genes = Vec::new();
        for r in rows {
            let f: Vec<&str> = r.split('\t').collect();
            assert_eq!(lof(f[g], f[hi], &hand), Some(f[lf]), "{}", f[g]);
            if f[g] == TTN {
                assert!(f[lim].starts_with("truncating_in chr2:") && f[lim].ends_with(TTN_NAME));
            }
            genes.push(f[g]);
        }
        assert_eq!(genes.len(), 84);
        for gene in hand.keys() {
            assert!(genes.contains(&gene.as_str()), "{gene}");
        }
    }
}
