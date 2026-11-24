//! `og dev build acmg-lof`: the `clingen_hi` and `lof` columns of
//! `assets/acmg_sf.tsv`, built from ClinGen's gene curation list.
//!
//! `scripts/ci/build-acmg-lof.sh` fetches the pinned list (an archived
//! copy, so its bytes never change) and checks its sha256; this verb only
//! reads it. It rewrites the two cells of each gene row and leaves every
//! other byte of the table as it is, so the table is its own oracle:
//! `--check` exits 1 when a rebuild would change it.
//!
//! `clingen_hi` is ClinGen's haploinsufficiency score as ClinGen writes it
//! (`3` sufficient evidence, `2` some, `1` little, `0` none, `30` a gene of
//! an autosomal recessive disease, `40` dosage sensitivity unlikely), or `.`
//! for a gene ClinGen has not curated. `lof` says whether a predicted
//! loss-of-function call (a VEP HIGH-impact record) counts as a secondary
//! finding in the gene (`og summary collect`, `collect/acmg.rs`):
//!
//! - `yes` for a dominant, semidominant or X-linked gene with score `3`, and
//!   for a recessive gene unless ClinGen scores it `40`: a recessive disease
//!   is the loss of both copies;
//! - `no` otherwise, so a PCSK9 stop-gain (score `40`; PCSK9 causes familial
//!   hypercholesterolemia by gain of function) is no finding;
//! - the hand rows below, each with its reason.

use super::{Args, BuildError, input_error};
use indexmap::IndexMap;
use std::io::Write;
use std::path::Path;

const USAGE: &str = "usage: og dev build acmg-lof --clingen FILE --table FILE [--check]\n";

/// Genes whose ClinGen score belongs to another disease than the one ACMG
/// lists them for. APOB's score 3 is for hypobetalipoproteinemia, while
/// ACMG lists APOB for familial hypercholesterolemia, which LDL receptor
/// binding variants cause. RET's score 3 is for Hirschsprung disease, while
/// ACMG lists RET for MEN2, a gain of function.
const NOT_LOF: [&str; 2] = ["APOB", "RET"];

/// Genes where a predicted loss-of-function call counts, with a mark:
/// Loeys-Dietz syndrome is mostly caused by missense variants of TGFBR1 and
/// TGFBR2, and truncating variants are a less certain cause.
const MARK: [&str; 2] = ["TGFBR1", "TGFBR2"];

/// Genes where ACMG reports truncating variants only (Table 1 footnote f
/// for TTN), whatever their score: the `limit` column of the table says
/// which truncating variants count.
const TRUNCATING: [&str; 1] = ["TTN"];

/// The `lof` cell of a gene.
pub fn lof(gene: &str, inheritance: &str, hi: &str) -> &'static str {
    if MARK.contains(&gene) {
        "mark"
    } else if TRUNCATING.contains(&gene) {
        "yes"
    } else if NOT_LOF.contains(&gene) {
        "no"
    } else if inheritance == "AR" {
        if hi == "40" { "no" } else { "yes" }
    } else if hi == "3" {
        "yes"
    } else {
        "no"
    }
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
            BuildError::Input(format!("{name}: a gene row comes before the #Gene Symbol line"))
        })?;
        let f: Vec<&str> = line.split('\t').collect();
        let (Some(gene), Some(score)) = (f.get(g), f.get(s)) else {
            return Err(BuildError::Input(format!("{name}: a row has too few columns: {line:?}")));
        };
        let score = score.trim();
        if gene.is_empty() {
            return Err(BuildError::Input(format!("{name}: a row has no gene symbol")));
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

/// The table with its `clingen_hi` and `lof` cells rebuilt from `scores`,
/// and the genes whose cells changed.
pub fn rebuild(
    table: &str,
    scores: &IndexMap<String, String>,
) -> Result<(String, Vec<String>), BuildError> {
    let mut out = String::with_capacity(table.len());
    let mut header: Option<[usize; 4]> = None;
    let mut changed = Vec::new();
    let mut seen = Vec::new();
    for line in table.split_inclusive('\n') {
        if line.starts_with('#') {
            out.push_str(line);
            continue;
        }
        let body = line.trim_end_matches(['\n', '\r']);
        let end = &line[body.len()..];
        let mut cells: Vec<String> = body.split('\t').map(str::to_owned).collect();
        let Some([g, inh, hi, lf]) = header else {
            let at = |c: &str| {
                cells.iter().position(|x| x == c).ok_or_else(|| {
                    BuildError::Input(format!("acmg_sf.tsv has no {c} column"))
                })
            };
            header = Some([at("gene")?, at("inheritance")?, at("clingen_hi")?, at("lof")?]);
            out.push_str(line);
            continue;
        };
        if body.is_empty() {
            out.push_str(line);
            continue;
        }
        if cells.len() <= [g, inh, hi, lf].into_iter().max().unwrap_or(0) {
            return Err(BuildError::Input(format!(
                "acmg_sf.tsv: a row has too few columns: {body:?}"
            )));
        }
        let gene = cells[g].clone();
        let score = scores.get(&gene).map_or(".", String::as_str).to_owned();
        let l = lof(&gene, &cells[inh], &score).to_owned();
        if cells[hi] != score || cells[lf] != l {
            changed.push(gene.clone());
        }
        cells[hi] = score;
        cells[lf] = l;
        seen.push(gene);
        out.push_str(&cells.join("\t"));
        out.push_str(end);
    }
    if header.is_none() {
        return Err(BuildError::Input("acmg_sf.tsv has no header row".into()));
    }
    for g in NOT_LOF.iter().chain(&MARK).chain(&TRUNCATING) {
        if !seen.iter().any(|s| s == g) {
            return Err(BuildError::Data(format!(
                "{g} has a hand rule in og dev build acmg-lof but no row in acmg_sf.tsv"
            )));
        }
    }
    Ok((out, changed))
}

fn run(args: &[String], err: &mut dyn Write) -> Result<u8, Result<BuildError, String>> {
    let a = Args::parse(args, &["--clingen", "--table"], &["--check"]).map_err(Err)?;
    let (clingen, table) = (
        a.need("--clingen").map_err(Err)?,
        a.need("--table").map_err(Err)?,
    );
    let scores = clingen_scores(Path::new(clingen)).map_err(Ok)?;
    let path = Path::new(table);
    let old = std::fs::read_to_string(path).map_err(|e| Ok(input_error(path, &e)))?;
    let (new, changed) = rebuild(&old, &scores).map_err(Ok)?;
    if a.has("--check") {
        if new == old {
            return Ok(0);
        }
        let _ = writeln!(
            err,
            "{table}: not what ClinGen's list gives: {} (run scripts/ci/build-acmg-lof.sh)",
            if changed.is_empty() {
                "the rows' bytes differ".to_owned()
            } else {
                format!("clingen_hi or lof differs for {}", changed.join(", "))
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
RET\t5979\t10q11.21\tchr10:1-2\t3\tSufficient evidence\n\
TGFBR1\t7046\t9q22.33\tchr9:1-2\t3\tSufficient evidence\n\
TGFBR2\t7048\t3p24.1\tchr3:1-2\t0\tNo evidence available\n\
TTN\t7273\t2q31.2\tchr2:3-4\t2\tSome evidence\n";

    const TABLE: &str = "# a comment\ngene\tinheritance\tnote\tclingen_hi\tlof\tlimit\n\
BRCA2\tAD\t.\t\t\t.\nPCSK9\tAD\t.\t\t\t.\nMUTYH\tAR\t.\t\t\t.\nMYH7\tAD\t.\t\t\t.\n\
APOB\tAD\t.\t\t\t.\nRET\tAD\t.\t\t\t.\nTGFBR1\tAD\t.\t\t\t.\nTGFBR2\tAD\t.\t\t\t.\n\
TTN\tAD\t.\t\t\t.\nGAA\tAR\t.\t\t\t.\nDES\tAD\t.\t\t\t.\n";

    fn scores() -> IndexMap<String, String> {
        let d = std::env::temp_dir().join(format!("og-acmg-lof-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let p = d.join("clingen.tsv");
        std::fs::write(&p, CLINGEN).unwrap();
        let s = clingen_scores(&p).unwrap();
        std::fs::remove_dir_all(&d).unwrap();
        s
    }

    fn cells(table: &str, gene: &str) -> (String, String) {
        let row = table
            .lines()
            .find(|l| l.starts_with(&format!("{gene}\t")))
            .unwrap();
        let f: Vec<&str> = row.split('\t').collect();
        (f[3].to_owned(), f[4].to_owned())
    }

    #[test]
    fn each_gene_gets_its_score_and_whether_a_loss_of_function_call_counts() {
        let (t, changed) = rebuild(TABLE, &scores()).unwrap();
        let want = [
            ("BRCA2", "3", "yes"),
            ("PCSK9", "40", "no"),
            ("MUTYH", "30", "yes"),
            ("MYH7", "0", "no"),
            ("APOB", "3", "no"),
            ("RET", "3", "no"),
            ("TGFBR1", "3", "mark"),
            ("TGFBR2", "0", "mark"),
            ("TTN", "2", "yes"),
            ("GAA", ".", "yes"),
            ("DES", ".", "no"),
        ];
        for (g, hi, l) in want {
            assert_eq!(cells(&t, g), (hi.to_owned(), l.to_owned()), "{g}");
        }
        assert_eq!(changed.len(), want.len());
        // The comment, the header and the other cells keep their bytes.
        assert!(t.starts_with("# a comment\ngene\tinheritance\tnote\tclingen_hi\tlof\tlimit\n"));
        assert!(t.contains("\nBRCA2\tAD\t.\t3\tyes\t.\n"));
    }

    #[test]
    fn a_rebuild_is_byte_stable_and_check_sees_a_changed_cell() {
        let s = scores();
        let (once, _) = rebuild(TABLE, &s).unwrap();
        let (twice, changed) = rebuild(&once, &s).unwrap();
        assert_eq!(once, twice);
        assert!(changed.is_empty());
        // CRLF endings stay CRLF.
        let crlf = once.replace('\n', "\r\n");
        assert_eq!(rebuild(&crlf, &s).unwrap().0, crlf);
        let edited = once.replace("PCSK9\tAD\t.\t40\tno", "PCSK9\tAD\t.\t40\tyes");
        assert_eq!(rebuild(&edited, &s).unwrap().1, ["PCSK9"]);

        let d = std::env::temp_dir().join(format!("og-acmg-lof-cli-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(d.join("clingen.tsv"), CLINGEN).unwrap();
        let table = d.join("acmg_sf.tsv");
        let args = |check: bool| {
            let mut a: Vec<String> = ["--clingen", "--table"]
                .iter()
                .zip([d.join("clingen.tsv"), table.clone()])
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
        assert!(e.contains("PCSK9"), "{e}");
        assert_eq!(call(false).0, 0);
        assert_eq!(std::fs::read_to_string(&table).unwrap(), once);
        assert_eq!(call(true), (0, String::new()));
        std::fs::write(&table, &edited).unwrap();
        let (code, e) = call(true);
        assert_eq!(code, 1);
        assert!(e.contains("lof differs for PCSK9"), "{e}");
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn bad_inputs_are_refused() {
        let s = scores();
        assert!(matches!(
            rebuild("gene\tinheritance\tlof\nA\tAD\tyes\n", &s),
            Err(BuildError::Input(_))
        ));
        // A hand rule whose gene left the table is stale.
        let no_ttn: String = TABLE
            .lines()
            .filter(|l| !l.starts_with("TTN"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert!(matches!(rebuild(&no_ttn, &s), Err(BuildError::Data(m)) if m.contains("TTN")));
        let d = std::env::temp_dir().join(format!("og-acmg-lof-bad-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
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
}
