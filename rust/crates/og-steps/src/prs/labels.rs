//! The score list, `assets/pgs_scores.tsv`: per PGS id its trait label, the
//! sex it applies to and its PGS Catalog evaluations.
//!
//! The file is tab separated with a header row (`pgs_id trait_reported
//! applies_to ppm_id ancestry effect`, the last four optional); lines that
//! start with `#` are comments. It is read as Python's `csv.DictReader`
//! reads it, so quoting and short rows behave as they did.

use og_io::{Lines, Utf8};
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::text::{is_decimal, lower, strip};
use std::collections::HashMap;
use std::io;
use std::path::Path;

/// The sexes a score can apply to.
pub const SEXES: [&str; 2] = ["female", "male"];
/// The 1000 Genomes groups pgsc_calc places a sample in.
pub const GROUPS: [&str; 5] = ["AFR", "AMR", "EAS", "EUR", "SAS"];

/// One PGS Catalog evaluation of a score.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Eval {
    /// One of the five 1000 Genomes groups, `AFR` to `SAS`.
    pub ancestry: String,
    /// `OR <number>` or `HR <number>`, as written.
    pub effect: String,
    /// `PPM` and six digits.
    pub ppm: String,
}

/// What the score list says about each PGS id.
#[derive(Debug, Default)]
pub struct Labels {
    traits: HashMap<String, String>,
    applies: HashMap<String, String>,
    evals: HashMap<String, Vec<Eval>>,
}

impl Labels {
    /// The trait label, when the list names one (an empty one counts as none).
    pub fn trait_of(&self, pid: &str) -> Option<&str> {
        self.traits
            .get(pid)
            .map(String::as_str)
            .filter(|t| !t.is_empty())
    }

    /// `female` or `male` for a score of one sex; `None` for both.
    pub fn applies_to(&self, pid: &str) -> Option<&str> {
        self.applies
            .get(pid)
            .map(String::as_str)
            .filter(|s| !s.is_empty())
    }

    /// The evaluations, in the list's order; empty when it names none.
    pub fn evals(&self, pid: &str) -> &[Eval] {
        self.evals.get(pid).map_or(&[], Vec::as_slice)
    }

    /// Reads the list at `path`. The error is the message after `ERROR: `.
    pub fn read(path: &Path) -> Result<Labels, String> {
        let name = path.file_name().map_or_else(
            || path.display().to_string(),
            |n| n.to_string_lossy().into_owned(),
        );
        let bytes =
            std::fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        if let Err(e) = std::str::from_utf8(&bytes) {
            return Err(format!("{name}: not UTF-8 at byte {}", e.valid_up_to()));
        }
        Labels::parse(&bytes, &name)
    }

    /// Parses the bytes of a list called `name` (UTF-8 already checked).
    fn parse(bytes: &[u8], name: &str) -> Result<Labels, String> {
        let mut lines = Lines::new(bytes, Utf8::Strict);
        let mut kept: Vec<io::Result<String>> = Vec::new();
        let mut line = String::new();
        loop {
            match lines.read_line(&mut line) {
                Ok(true) if line.starts_with('#') => {}
                Ok(true) => kept.push(Ok(std::mem::take(&mut line))),
                Ok(false) => break,
                Err(e) => return Err(format!("{name}: {e}")),
            }
        }
        let rows = DictReader::new(kept.into_iter(), '\t', Quoting::Minimal)
            .collect::<Result<Vec<DictRow>, _>>()
            .map_err(|e| format!("{name}: {e}"))?;
        // The checks run in the order the Python ran them: every row's
        // applies_to, then every row's evaluations.
        let mut out = Labels::default();
        for r in &rows {
            let pid = strip(cell(r, "pgs_id"), None);
            let sex = lower(strip(cell(r, "applies_to"), None));
            if pid.is_empty() {
                continue;
            }
            if !sex.is_empty() && !SEXES.contains(&sex.as_str()) {
                return Err(format!(
                    "{name}: applies_to of {pid} is '{sex}', expected female, male or empty"
                ));
            }
            out.applies.insert(pid.to_owned(), sex);
        }
        for r in &rows {
            read_evals(r, name, &mut out.evals)?;
        }
        for r in &rows {
            let pid = cell(r, "pgs_id");
            if !pid.is_empty() {
                let label = strip(cell(r, "trait_reported"), None);
                out.traits
                    .insert(strip(pid, None).to_owned(), label.to_owned());
            }
        }
        Ok(out)
    }
}

/// `r.get(key) or ""`.
fn cell<'a>(r: &'a DictRow, key: &str) -> &'a str {
    r.get(key).flatten().unwrap_or("")
}

/// A `;` list of a row: its stripped entries, or none for a blank cell.
fn entries<'a>(r: &'a DictRow, key: &str) -> Vec<&'a str> {
    let v = cell(r, key);
    if strip(v, None).is_empty() {
        Vec::new()
    } else {
        v.split(';').map(|x| strip(x, None)).collect()
    }
}

/// The evaluations of one row into `out`, refused as the Python refused them.
fn read_evals(r: &DictRow, name: &str, out: &mut HashMap<String, Vec<Eval>>) -> Result<(), String> {
    let pid = strip(cell(r, "pgs_id"), None);
    let (ppm, anc, eff) = (
        entries(r, "ppm_id"),
        entries(r, "ancestry"),
        entries(r, "effect"),
    );
    if ppm.is_empty() && anc.is_empty() && eff.is_empty() {
        return Ok(());
    }
    if pid.is_empty() {
        return Err(format!(
            "{name}: a row has ppm_id, ancestry or effect entries but no pgs_id"
        ));
    }
    if out.contains_key(pid) {
        return Err(format!("{name}: {pid} has evaluations on two rows"));
    }
    if !(ppm.len() == anc.len() && anc.len() == eff.len()) {
        return Err(format!(
            "{name}: {pid} has {} ppm_id, {} ancestry and {} effect entries",
            ppm.len(),
            anc.len(),
            eff.len()
        ));
    }
    let mut evals = Vec::with_capacity(ppm.len());
    for ((p, a), e) in ppm.iter().zip(&anc).zip(&eff) {
        if !GROUPS.contains(a) {
            return Err(format!(
                "{name}: ancestry of {pid} is '{a}', expected one of {}",
                GROUPS.join(", ")
            ));
        }
        if !is_ppm_id(p) {
            return Err(format!(
                "{name}: ppm_id of {pid} is '{p}', expected PPM and six digits"
            ));
        }
        if !is_effect(e) {
            return Err(format!(
                "{name}: effect of {pid} is '{e}', expected 'OR <number>' or 'HR <number>'"
            ));
        }
        evals.push(Eval {
            ancestry: (*a).to_owned(),
            effect: (*e).to_owned(),
            ppm: (*p).to_owned(),
        });
    }
    for (i, e) in evals.iter().enumerate() {
        if evals[..i].iter().any(|x| x.ancestry == e.ancestry) {
            return Err(format!("{name}: {pid} lists an ancestry twice"));
        }
    }
    out.insert(pid.to_owned(), evals);
    Ok(())
}

/// `\d+` of Python's `re` on a `str`: one or more Unicode decimal digits.
fn digits(s: &str) -> bool {
    !s.is_empty() && s.chars().all(is_decimal)
}

/// `re.fullmatch(r"PPM\d{6}", s)`.
fn is_ppm_id(s: &str) -> bool {
    s.strip_prefix("PPM")
        .is_some_and(|d| d.chars().count() == 6 && digits(d))
}

/// `re.fullmatch(r"(OR|HR) \d+(\.\d+)?", s)`.
fn is_effect(s: &str) -> bool {
    let Some(num) = s.strip_prefix("OR ").or_else(|| s.strip_prefix("HR ")) else {
        return false;
    };
    match num.split_once('.') {
        None => digits(num),
        Some((int, frac)) => digits(int) && digits(frac),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &str) -> Result<Labels, String> {
        Labels::parse(text.as_bytes(), "list.tsv")
    }

    const HEAD: &str = "pgs_id\ttrait_reported\tapplies_to\tppm_id\tancestry\teffect\n";

    #[test]
    fn reads_traits_sexes_and_evaluations() {
        let l = parse(&format!(
            "# a comment\n{HEAD}PGS1\t Coronary artery disease \tFemale\tPPM000001;PPM000002\tEUR;AFR\tHR 1.53;OR 2\n\
             PGS2\t\t\t\t\t\n"
        ))
        .expect("read");
        assert_eq!(l.trait_of("PGS1"), Some("Coronary artery disease"));
        assert_eq!(l.trait_of("PGS2"), None);
        assert_eq!(l.applies_to("PGS1"), Some("female"));
        assert_eq!(l.applies_to("PGS2"), None);
        let e = l.evals("PGS1");
        assert_eq!(e.len(), 2);
        assert_eq!(
            (
                e[1].ancestry.as_str(),
                e[1].effect.as_str(),
                e[1].ppm.as_str()
            ),
            ("AFR", "OR 2", "PPM000002")
        );
        assert!(l.evals("PGS2").is_empty());
    }

    #[test]
    fn refuses_what_the_python_refused() {
        for (row, want) in [
            ("PGS1\tx\tboth\t\t\t", "applies_to of PGS1 is 'both'"),
            (
                "PGS1\tx\t\tPPM000001;PPM000002\tEUR\tHR 1",
                "has 2 ppm_id, 1 ancestry and 1 effect entries",
            ),
            (
                "PGS1\tx\t\tPPM000001\tEuropean\tHR 1",
                "ancestry of PGS1 is 'European'",
            ),
            (
                "PGS1\tx\t\tPPM000001;PPM000002\tEUR;EUR\tHR 1;HR 2",
                "lists an ancestry twice",
            ),
            (
                "PGS1\tx\t\tPGS000001\tEUR\tHR 1",
                "ppm_id of PGS1 is 'PGS000001'",
            ),
            (
                "PGS1\tx\t\tPPM000001\tEUR\tbeta 1",
                "effect of PGS1 is 'beta 1'",
            ),
            (
                "PGS1\tx\t\tPPM000001\tEUR\tOR 1.",
                "effect of PGS1 is 'OR 1.'",
            ),
            ("\tx\t\tPPM000001\tEUR\tHR 1", "but no pgs_id"),
            (
                "PGS1\tx\t\tPPM000001\tEUR\tHR 1\nPGS1\tx\t\tPPM000002\tAFR\tHR 2",
                "evaluations on two rows",
            ),
        ] {
            let got = parse(&format!("{HEAD}{row}\n")).err().unwrap_or_default();
            assert!(
                got.starts_with("list.tsv: ") && got.contains(want),
                "{row}: {got}"
            );
        }
    }

    #[test]
    fn a_bad_sex_is_reported_before_a_bad_evaluation_of_an_earlier_row() {
        let got = parse(&format!(
            "{HEAD}PGS1\tx\t\tPPM1\tEUR\tHR 1\nPGS2\tx\tnone\t\t\t\n"
        ));
        assert!(got.err().unwrap_or_default().contains("applies_to of PGS2"));
    }

    #[test]
    fn unicode_digits_count_as_python_counts_them() {
        assert!(is_ppm_id("PPM00000\u{0663}"));
        assert!(!is_ppm_id("PPM00000"));
        assert!(is_effect("HR \u{0661}.5"));
        assert!(!is_effect("HR 1.5.1"));
    }

    #[test]
    fn a_file_that_is_not_utf8_is_a_message() {
        let dir = std::env::temp_dir().join(format!("og-prs-labels-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let p = dir.join("latin1.tsv");
        std::fs::write(&p, b"pgs_id\ttrait_reported\nPGS1\tcaf\xe9\n").expect("write");
        assert_eq!(
            Labels::read(&p).err().as_deref(),
            Some("latin1.tsv: not UTF-8 at byte 30")
        );
        std::fs::remove_dir_all(&dir).expect("clean");
    }
}
