//! The eight small installers og runs: `--pull-only`, `--sample-qc-data`,
//! `--viral-data`, `--kir-data`, `--chip-data`, `--cnvpytor-data`,
//! `--somatic-data` and `--prs-academic`. Each does what its bash function
//! in scripts/setup.sh does, with the same downloads (through
//! [`Run::fetch`], og setup fetch in process), the same files, the same
//! messages and the same exit status; tests/setup/ab.sh compares the trees
//! they leave with bash's.

use super::{Run, nonempty, sha256};
use crate::images::{self, Engine};
use crate::lock;
use std::io::Write;
use std::path::{Path, PathBuf};

/// The keys this module runs.
pub const PORTED: &[&str] = &[
    "pull-only",
    "sample-qc-data",
    "viral-data",
    "kir-data",
    "chip-data",
    "cnvpytor-data",
    "somatic-data",
    "prs-academic",
];

/// somalier's sites for step 33; setup.sh holds the same pin while it has
/// its bash installer (tests/core.rs compares them).
pub const SOMALIER_SITES_URL: &str =
    "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz";
pub const SOMALIER_SITES_SHA256: &str =
    "d1a853b8bb2e5f1a520bc67c8303be699543d225d785bce51e8335a2420489b5";
/// The VerifyBamID2 panel of the release VERIFYBAMID2_IMAGE packages
/// (v2.0.3): each file's URL and sha256.
pub const VB2_PANEL_FILES: [(&str, &str); 3] = [
    (
        "https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD",
        "259e320123756bb702542e3b3c4d766b481426d96cd63abd4b9f9fced035f375",
    ),
    (
        "https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu",
        "f7e8b4fad17cc433d887f18bf183fa4b89597d5a375f22133f44f39c5e35ffc8",
    ),
    (
        "https://raw.githubusercontent.com/Griffan/VerifyBamID/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed",
        "0025d782137e5906bbc7afc553c49d2c82c5b8f0bfa0649bc415887442f42869",
    ),
];

/// Where the PGS Catalog keeps its scores: PGS_BASE_URL, or this.
pub const PGS_BASE_URL: &str = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores";

/// Runs KEY when it is one of [`PORTED`]; its exit status.
pub fn run(key: &str, r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    Some(match key {
        "pull-only" => pull_only(r, out, err),
        "sample-qc-data" => sample_qc_data(r, out, err),
        "viral-data" => viral_data(r, out, err),
        "kir-data" => kir_data(r, out, err),
        "chip-data" => chip_data(r, out, err),
        "cnvpytor-data" => cnvpytor_data(r, out, err),
        "somatic-data" => somatic_data(r, out, err),
        "prs-academic" => prs_academic(r, out, err),
        _ => return None,
    })
}

fn p(s: &str) -> PathBuf {
    PathBuf::from(s)
}

fn show(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// `setup.sh --pull-only`: the header, the engine check, and every image
/// of versions.env that is not marked optional, pulled and tagged.
pub fn pull_only(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let g = if r.genome.is_empty() {
        "(none: --pull-only)"
    } else {
        r.genome.as_str()
    };
    let _ = write!(
        out,
        "============================================\n  OwnGenome — Setup\n  Data directory: {g}\n============================================\n\n"
    );
    let e = Engine::from_vars(&r.vars);
    if !super::command_v(&e.cmd) {
        let _ = writeln!(out, "ERROR: Docker is not installed.");
        let _ = writeln!(out, "  Install Docker: https://docs.docker.com/get-docker/");
        return 1;
    }
    if !e.quiet(&["info"]) {
        let _ = writeln!(out, "ERROR: Docker daemon is not running.");
        let _ = writeln!(
            out,
            "  Start Docker Desktop or run: sudo systemctl start docker"
        );
        return 1;
    }
    let _ = writeln!(out, "[OK] Docker is running.");
    let _ = writeln!(out);
    let _ = writeln!(out, "=== Docker Images (~10-15 GB total) ===");
    let list = match images::pipeline_images(&r.root) {
        Ok(l) => l,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            return 1;
        }
    };
    u8::from(!images::pull_images(&e, &list, &r.vars, out, err))
}

/// install_sample_qc_data: somalier's sites and VerifyBamID2's panel under
/// reference/, each checked against its sha256 on every run; a file there
/// with another sha256 is downloaded again.
pub fn sample_qc_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let dest = format!("{}/reference", r.genome);
    let mut rc = 0;
    let f = p(&format!("{dest}/somalier/sites.hg38.vcf.gz"));
    if nonempty(&f) && sha256(&f) == SOMALIER_SITES_SHA256 {
        let _ = writeln!(
            out,
            "[OK] somalier sites (step 33) already present (sha256 matches)."
        );
    } else {
        if nonempty(&f) {
            let _ = writeln!(
                out,
                "[WARN] {} is not the pinned file (sha256 differs): downloading it again.",
                show(&f)
            );
        }
        if r.fetch(SOMALIER_SITES_URL, &f, "sha256", SOMALIER_SITES_SHA256, err) {
            let _ = writeln!(out, "[OK] somalier sites (step 33): {}", show(&f));
        } else {
            let _ = writeln!(
                out,
                "[WARN] Could not install somalier's sites. Step 33 needs them; re-run setup.sh."
            );
            rc = 1;
        }
    }
    for (url, sha) in VB2_PANEL_FILES {
        let name = url.rsplit('/').next().unwrap_or(url);
        let f = p(&format!("{dest}/verifybamid2/{name}"));
        if nonempty(&f) && sha256(&f) == sha {
            continue;
        }
        if nonempty(&f) {
            let _ = writeln!(
                out,
                "[WARN] {} is not the pinned file (sha256 differs): downloading it again.",
                show(&f)
            );
        }
        if !r.fetch(url, &f, "sha256", sha, err) {
            let _ = writeln!(
                out,
                "[WARN] Could not install {name} of VerifyBamID2's panel. Step 33 needs it; re-run setup.sh."
            );
            rc = 1;
        }
    }
    if rc == 0 {
        let _ = writeln!(
            out,
            "[OK] VerifyBamID2 panel (step 33): {dest}/verifybamid2/"
        );
    }
    rc
}

/// install_viral_data: the viral genomes of step 47 (the table's
/// viral_genomes row), checked against VIRAL_REF_SHA256.
pub fn viral_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let v = r.version("VIRAL_REF_VERSION");
    let dest = p(&format!(
        "{}/reference/viral-{v}/viral_genomes.fa",
        r.genome
    ));
    let want = r.version("VIRAL_REF_SHA256");
    if nonempty(&dest) && sha256(&dest) == want {
        let _ = writeln!(out, "[OK] Viral genomes {v} (step 47) already present.");
        return 0;
    }
    // Not kept: an installed file that fails its sha256.
    let _ = std::fs::remove_file(&dest);
    let Some(url) = r.source_of("viral_genomes", "url", err) else {
        return 1;
    };
    if !r.fetch(&url, &dest, "sha256", &want, err) {
        let _ = writeln!(
            out,
            "[WARN] Could not install the viral genomes of step 47 (NCBI E-utilities)."
        );
        return 1;
    }
    let n = std::fs::read(&dest)
        .map(|b| {
            b.split(|&c| c == b'\n')
                .filter(|l| l.first() == Some(&b'>'))
                .count()
        })
        .unwrap_or(0);
    let _ = writeln!(
        out,
        "[OK] Viral genomes {v} (step 47, {n} sequences): {}",
        show(&dest)
    );
    0
}

/// Whether a line of TEXT matches PATTERN as `grep -q` reads a pattern of
/// letters, digits, spaces, '-' and '.' (a '.' matches any one character).
pub fn grep_line(text: &[u8], pattern: &str) -> bool {
    let pat = pattern.as_bytes();
    text.split(|&c| c == b'\n').any(|line| {
        line.len() >= pat.len()
            && (0..=line.len() - pat.len()).any(|i| {
                pat.iter()
                    .zip(&line[i..i + pat.len()])
                    .all(|(&p, &c)| p == b'.' || p == c)
            })
    })
}

/// The md5 of kir.dat in a release's md5checksum.txt ("MD5 (kir.dat) =
/// <md5>"): the last word of each line whose second word is "(kir.dat)".
pub fn kir_md5(list: &[u8]) -> String {
    String::from_utf8_lossy(list)
        .lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split([' ', '\t']).filter(|w| !w.is_empty()).collect();
            (w.get(1) == Some(&"(kir.dat)")).then(|| w.last().copied().unwrap_or("").to_owned())
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// install_kir_data: IPD-KIR KIR_DB_RELEASE's kir.dat (step 08, KIR=true),
/// checked against the md5 its release lists, and refused when it names
/// another release.
pub fn kir_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let rel = r.version("KIR_DB_RELEASE");
    let branch = rel.replace('.', "");
    let dest = p(&format!("{}/kir/IPD-KIR_{rel}/kir.dat", r.genome));
    if nonempty(&dest) {
        let _ = writeln!(
            out,
            "[OK] IPD-KIR {rel} (step 08, KIR=true) already present."
        );
        return 0;
    }
    let url = format!("https://raw.githubusercontent.com/ANHIG/IPDKIR/{branch}");
    let md5 = r
        .get_quiet(&format!("{url}/md5checksum.txt"))
        .map(|b| kir_md5(&b))
        .unwrap_or_default();
    if md5.is_empty() || !r.fetch(&format!("{url}/kir.dat"), &dest, "md5", &md5, err) {
        let _ = writeln!(
            out,
            "[WARN] Could not install IPD-KIR {rel} (is it a release with a branch {branch}?)."
        );
        return 1;
    }
    let says = std::fs::read(&dest)
        .is_ok_and(|b| grep_line(&b, &format!("IPD-KIR Release Version {rel}")));
    if !says {
        // Not kept: the wrong release, which step 08 must not read.
        let _ = std::fs::remove_file(&dest);
        let _ = writeln!(
            out,
            "[WARN] {url}/kir.dat does not say IPD-KIR Release Version {rel}; removed it."
        );
        return 1;
    }
    let _ = writeln!(
        out,
        "[OK] IPD-KIR {rel} (step 08, KIR=true): {}",
        show(&dest)
    );
    0
}

/// install_chip_data: Ensembl ENSEMBL_GFF_RELEASE's GRCh38 GFF3 (step 49),
/// checked against ENSEMBL_GFF_SHA256.
pub fn chip_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let Some(dest) = r.source_of("ensembl_gff", "path", err) else {
        return 1;
    };
    let Some(url) = r.source_of("ensembl_gff", "url", err) else {
        return 1;
    };
    let dest = p(&dest);
    let rel = r.version("ENSEMBL_GFF_RELEASE");
    let want = r.version("ENSEMBL_GFF_SHA256");
    if nonempty(&dest) && sha256(&dest) == want {
        let _ = writeln!(
            out,
            "[OK] Ensembl {rel} gene annotation (step 49) already present."
        );
        return 0;
    }
    // Not kept: an installed file that fails its sha256.
    let _ = std::fs::remove_file(&dest);
    if !r.fetch(&url, &dest, "sha256", &want, err) {
        let _ = writeln!(
            out,
            "[WARN] Could not download Ensembl {rel}'s gene annotation. Step 49 needs it; run setup.sh --chip-data again."
        );
        return 1;
    }
    let _ = writeln!(
        out,
        "[OK] Ensembl {rel} gene annotation (step 49): {}",
        show(&dest)
    );
    0
}

/// setup_pin_rows: the versions.env value of each pin name KEYS lists, one
/// row each ("FILE KIND VALUE").
pub fn pin_rows(r: &Run, keys: &str) -> Vec<String> {
    keys.split_whitespace().map(|k| r.version(k)).collect()
}

/// `read -r a b c` of a row: the first two words and the rest.
pub fn three(row: &str) -> (String, String, String) {
    let t = row.trim_matches([' ', '\t']);
    let (a, rest) = t.split_once([' ', '\t']).unwrap_or((t, ""));
    let rest = rest.trim_start_matches([' ', '\t']);
    let (b, c) = rest.split_once([' ', '\t']).unwrap_or((rest, ""));
    (
        a.to_owned(),
        b.to_owned(),
        c.trim_matches([' ', '\t']).to_owned(),
    )
}

/// `mkdir -p DIR`, with mkdir's message when it fails.
fn mkdir(dir: &Path, err: &mut dyn Write) -> bool {
    match std::fs::create_dir_all(dir) {
        Ok(()) => true,
        Err(e) => {
            let _ = writeln!(err, "mkdir: cannot create directory '{}': {e}", show(dir));
            false
        }
    }
}

fn hold(lock_path: &Path, out: &mut dyn Write, err: &mut dyn Write) -> Option<lock::Lock> {
    match lock::acquire(lock_path, out) {
        Ok(l) => Some(l),
        Err(e) => {
            let _ = writeln!(err, "ERROR: cannot write {}/pid: {e}", show(lock_path));
            None
        }
    }
}

/// install_cnvpytor_data: CNVpytor's reference files (CNVPYTOR_DATA_FILES)
/// under reference/cnvpytor, under its lock; a file there is kept.
pub fn cnvpytor_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let rows = pin_rows(r, &r.version("CNVPYTOR_DATA_FILES"));
    let base = r.version("CNVPYTOR_DATA_BASE_URL");
    let dir = p(&format!("{}/reference/cnvpytor", r.genome));
    if !mkdir(&dir, err) {
        return 1;
    }
    let Some(_held) = hold(&dir.join(".setup.lock"), out, err) else {
        return 1;
    };
    let mut rc = 0;
    for row in rows_or_one(&rows) {
        let (file, kind, pin) = three(row);
        let dest = p(&format!("{}/{file}", show(&dir)));
        if nonempty(&dest) {
            continue;
        }
        if !r.fetch(&format!("{base}/{file}"), &dest, &kind, &pin, err) {
            rc = 1;
        }
    }
    rc
}

/// A here-string of no rows still gives `read` one empty line.
fn rows_or_one(rows: &[String]) -> Vec<&str> {
    if rows.is_empty() {
        vec![""]
    } else {
        rows.iter().map(String::as_str).collect()
    }
}

/// install_indexed_downloads DIR BASE KEYS: published data and index pairs
/// (rows two by two), under DIR's lock. A pair that is there is kept; an
/// index without its data is fetched again after the data; a failed index
/// leaves the data for a retry.
pub fn indexed_downloads(
    r: &Run,
    dir: &Path,
    base: &str,
    keys: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    let rows = pin_rows(r, keys);
    if !mkdir(dir, err) {
        return 1;
    }
    let Some(_held) = hold(&dir.join(".setup-indexed.lock"), out, err) else {
        return 1;
    };
    let rows = rows_or_one(&rows);
    let mut rc = 0;
    for pair in rows.chunks_exact(2) {
        let (file, kind, pin) = three(pair[0]);
        let (index, index_kind, index_pin) = three(pair[1]);
        let d = show(dir);
        let dest = p(&format!("{d}/{file}"));
        let idx = p(&format!("{d}/{index}"));
        if nonempty(&dest) && nonempty(&idx) {
            continue;
        }
        // Not kept: an index without its data, fetched again below.
        if let Err(e) = std::fs::remove_file(&idx)
            && e.kind() != std::io::ErrorKind::NotFound
        {
            let _ = writeln!(err, "rm: cannot remove '{}': {e}", show(&idx));
            return 1;
        }
        if !r.fetch(&format!("{base}/{file}"), &dest, &kind, &pin, err) {
            rc = 1;
            continue;
        }
        if !r.fetch(
            &format!("{base}/{index}"),
            &idx,
            &index_kind,
            &index_pin,
            err,
        ) {
            rc = 1;
        }
    }
    rc
}

/// install_somatic_data: GATK's somatic resources (SOMATIC_DATA_FILES)
/// under somatic/.
pub fn somatic_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    indexed_downloads(
        r,
        &p(&format!("{}/somatic", r.genome)),
        &r.version("SOMATIC_DATA_BASE_URL"),
        &r.version("SOMATIC_DATA_FILES"),
        out,
        err,
    )
}

/// The scores of LIST (assets/pgs_scores.tsv) whose licence column is WANT
/// (an empty cell is `author`), in the list's order.
pub fn pgs_ids(list: &str, want: &str) -> Vec<String> {
    let mut cols: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut header = false;
    let mut ids = Vec::new();
    for line in list.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let first = f.first().copied().unwrap_or("");
        if !header && first == "pgs_id" {
            for (i, name) in f.iter().enumerate() {
                cols.insert((*name).to_owned(), i);
            }
            header = true;
            continue;
        }
        let is_id = first
            .strip_prefix("PGS")
            .is_some_and(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()));
        if header && is_id {
            let l = cols
                .get("licence")
                .and_then(|&i| f.get(i).copied())
                .unwrap_or("");
            let l = if l.is_empty() { "author" } else { l };
            if l == want {
                ids.push(first.to_owned());
            }
        }
    }
    ids
}

/// fetch_prs_scores LICENCE: the GRCh38-harmonised file of each score of
/// assets/pgs_scores.tsv under LICENCE, each checked against the md5 the
/// PGS Catalog publishes beside it; a score there is kept.
pub fn prs_scores(r: &Run, licence: &str, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let list = r.root.join("assets/pgs_scores.tsv");
    let Ok(text) = std::fs::read_to_string(&list) else {
        let _ = writeln!(
            out,
            "[WARN] Could not read the PGS score list {}",
            show(&list)
        );
        return 1;
    };
    let dir = p(&format!("{}/prs_scores", r.genome));
    if !mkdir(&dir, err) {
        return 1;
    }
    let base = r
        .vars
        .nonempty("PGS_BASE_URL")
        .unwrap_or(PGS_BASE_URL)
        .to_owned();
    let mut rc = 0;
    for id in pgs_ids(&text, licence) {
        let dest = p(&format!("{}/prs_scores/{id}.txt.gz", r.genome));
        if nonempty(&dest) {
            continue;
        }
        let url = format!("{base}/{id}/ScoringFiles/Harmonized/{id}_hmPOS_GRCh38.txt.gz");
        if !r.fetch(&url, &dest, "md5", &format!("{url}.md5"), err) {
            let _ = writeln!(
                out,
                "[WARN] Could not download {id} (step 25 tries again when it runs): {url}"
            );
            rc = 1;
        }
    }
    rc
}

/// `--prs-academic`: the licence notices (not in an --everything row), then
/// the four scores free for academic research use only.
pub fn prs_academic(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if r.vars.get("SETUP_EVERYTHING_SCOPE") != Some("row") {
        for t in crate::components::PRS_ACADEMIC_TERMS {
            let _ = writeln!(out, "  {t}");
        }
    }
    if prs_scores(r, "academic-research", out, err) != 0 {
        return 1;
    }
    let _ = writeln!(
        out,
        "[OK] Academic-only PGS Catalog scores (step 25 with PRS_ACADEMIC=true): {}/prs_scores/",
        r.genome
    );
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_kir_md5_is_the_last_word_of_its_line() {
        let list = b"MD5 (hla.dat) = 111\nMD5 (kir.dat) = abc123\n";
        assert_eq!(kir_md5(list), "abc123");
        assert_eq!(kir_md5(b"MD5 (other) = 1\n"), "");
    }

    #[test]
    fn a_dot_of_the_release_matches_any_character() {
        assert!(grep_line(
            b"x\nIPD-KIR Release Version 2.15.0\n",
            "IPD-KIR Release Version 2.15.0"
        ));
        assert!(grep_line(
            b"IPD-KIR Release Version 2x15y0",
            "IPD-KIR Release Version 2.15.0"
        ));
        assert!(!grep_line(
            b"IPD-KIR Release Version 2.14.0",
            "IPD-KIR Release Version 2.15.0"
        ));
    }

    #[test]
    fn rows_read_as_read_reads_them() {
        assert_eq!(
            three("  a.vcf.gz  sha256 abc def "),
            ("a.vcf.gz".into(), "sha256".into(), "abc def".into())
        );
        assert_eq!(three(""), (String::new(), String::new(), String::new()));
    }

    #[test]
    fn the_academic_scores_are_picked_by_their_licence() {
        let list = "# comment\tx\npgs_id\tname\tlicence\nPGS000001\ta\t\nPGS000014\tb\tacademic-research\nPGSX\tc\tacademic-research\n";
        assert_eq!(pgs_ids(list, "academic-research"), vec!["PGS000014"]);
        assert_eq!(pgs_ids(list, "author"), vec!["PGS000001"]);
    }
}
