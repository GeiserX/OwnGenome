//! `og longread clair3-model`: the Clair3 model that matches the basecaller
//! of an Oxford Nanopore sample.
//!
//! Clair3 has one model per pore, basecaller accuracy and basecaller model
//! version. A model trained on other basecalls still runs and still writes a
//! VCF, but every call is less accurate, and nothing in the output says so.
//! So step 03e never assumes a basecaller: it reads the basecaller model
//! from the BAM's `@RG` `DS` field (`basecall_model=...`, which Dorado writes
//! and step 02b carries), or takes it from `BASECALL_MODEL`, maps it here to
//! a model directory, and stops when no model matches or the image lacks it.
//! `CLAIR3_MODEL` still names a model directory outright.
//!
//! The verb runs inside the Clair3 image, so it checks the directory on the
//! filesystem it sees, and lists the image's models when it stops.
//!
//! `--basecaller-only` prints the basecaller model of a BAM's header (step
//! 02b reads an unaligned BAM's this way) and maps nothing.

use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

use crate::paraphase::bam_header;

const USAGE: &str = "usage: og longread clair3-model [--bam BAM] [--basecaller NAME] [--model DIR] [--models-dir DIR] [--basecaller-only]
  Prints the Clair3 model directory for an ONT sample. --model (CLAIR3_MODEL)
  wins; else --basecaller (BASECALL_MODEL); else the basecall_model= of the
  BAM's @RG DS fields. Exit 1 when no model fits or the directory is missing.
  --basecaller-only prints the BAM's basecaller model (nothing when none).
";

/// The Clair3 model directory name for a basecaller model name, or `None`
/// when Clair3 has no model for it.
///
/// Dorado names a model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`; Clair3 names
/// the matching one `r1041_e82_400bps_sup_v500`, one digit per version part.
/// R9.4.1 basecalls (Dorado `dna_r9.4.1_e8_sup@v3.3`, Guppy
/// `dna_r9.4.1_450bps_sup`) take Clair3's two R9.4.1 models, trained on
/// Guppy 5 SUP and on Guppy 3.6/4.2 HAC. Fast basecalls, RNA, other pores
/// and other speeds have no Clair3 model.
pub fn clair3_model_for(basecaller: &str) -> Option<String> {
    let name = basecaller.trim();
    let name = name.strip_prefix("dna_").unwrap_or(name);
    let (config, version) = match name.split_once('@') {
        Some((c, v)) => (c, Some(v)),
        None => (name, None),
    };
    if let Some(rest) = config.strip_prefix("r10.4.1_e8.2_400bps_") {
        let accuracy = match rest {
            "sup" | "hac" => rest,
            _ => return None,
        };
        let digits = version_digits(version?)?;
        return Some(format!("r1041_e82_400bps_{accuracy}_v{digits}"));
    }
    if let Some(rest) = config.strip_prefix("r9.4.1_") {
        // e8_sup, 450bps_sup, 450bps_hac_prom: the pore's speed and the
        // Guppy "prom" suffix do not change which Clair3 model fits.
        let words: Vec<&str> = rest.split('_').collect();
        if words.contains(&"sup") {
            return Some("r941_prom_sup_g5014".to_owned());
        }
        if words.contains(&"hac") {
            return Some("r941_prom_hac_g360+g422".to_owned());
        }
    }
    None
}

/// `v5.0.0` -> `500`, `v4.1.0` -> `410`, `v5.0` -> `500`: one digit per
/// part, as Clair3's directory names have them. `None` for anything else.
fn version_digits(version: &str) -> Option<String> {
    let parts: Vec<&str> = version.strip_prefix('v')?.split('.').collect();
    if !(2..=3).contains(&parts.len()) {
        return None;
    }
    let mut digits = String::new();
    for p in &parts {
        if p.len() != 1 || !p.bytes().all(|b| b.is_ascii_digit()) {
            return None;
        }
        digits.push_str(p);
    }
    if parts.len() == 2 {
        digits.push('0');
    }
    Some(digits)
}

/// The distinct `basecall_model=` values of the header's `@RG` `DS`
/// fields, sorted, and how many read groups name none.
pub fn basecallers(header: &str) -> (Vec<String>, usize) {
    let mut found = BTreeSet::new();
    let mut unnamed = 0;
    for line in header.split('\n').filter(|l| l.starts_with("@RG\t")) {
        let mut named = false;
        for tag in line.trim_end_matches('\r').split('\t').skip(1) {
            let Some(ds) = tag.strip_prefix("DS:") else {
                continue;
            };
            for word in ds.split(' ') {
                if let Some(v) = word.strip_prefix("basecall_model=")
                    && !v.is_empty()
                {
                    found.insert(v.to_owned());
                    named = true;
                }
            }
        }
        if !named {
            unnamed += 1;
        }
    }
    (found.into_iter().collect(), unnamed)
}

/// The one basecaller model of every read group of the header; `None`
/// when no read group names one. Read groups that name different models,
/// or some that name one and some none, are an error: one Clair3 model
/// cannot be chosen for the whole BAM from them.
pub fn bam_basecaller(header: &str) -> Result<Option<String>, String> {
    let (found, unnamed) = basecallers(header);
    match found.as_slice() {
        [] => Ok(None),
        [one] if unnamed == 0 => Ok(Some(one.clone())),
        [one] => Err(format!(
            "{unnamed} read group(s) name no basecaller model while others name {one}, so the model of every read is not known"
        )),
        many => Err(format!(
            "the read groups name more than one basecaller model ({}), and one Clair3 model cannot fit them all",
            many.join(", ")
        )),
    }
}

/// The directories under `dir`, sorted; empty when it cannot be read.
fn models_in(dir: &str) -> Vec<String> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut names: Vec<String> = entries
        .flatten()
        .filter(|e| e.path().is_dir())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[derive(Default)]
struct Opts {
    bam: Option<String>,
    basecaller: Option<String>,
    model: Option<String>,
    models_dir: Option<String>,
    basecaller_only: bool,
}

fn parse(args: &[String]) -> Result<Opts, String> {
    let mut o = Opts::default();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let (key, inline) = match a.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let slot = match key {
            "--basecaller-only" if inline.is_none() => {
                o.basecaller_only = true;
                continue;
            }
            "--bam" => &mut o.bam,
            "--basecaller" => &mut o.basecaller,
            "--model" => &mut o.model,
            "--models-dir" => &mut o.models_dir,
            _ => return Err(format!("unrecognized arguments: {a}")),
        };
        let value = match inline {
            Some(v) => v,
            None => it
                .next()
                .ok_or_else(|| format!("argument {key}: expected one argument"))?
                .clone(),
        };
        *slot = Some(value);
    }
    Ok(o)
}

/// `og longread clair3-model ARGS`: the model directory on stdout, or the
/// reason on stderr and exit 1 (2 for a usage error).
pub fn clair3_model_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = out.write_all(USAGE.as_bytes());
        return 0;
    }
    let o = match parse(args) {
        Ok(o) => o,
        Err(msg) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og longread clair3-model: error: {msg}");
            return 2;
        }
    };
    match run(&o) {
        Ok(Some(line)) => {
            let _ = writeln!(out, "{line}");
            0
        }
        Ok(None) => 0,
        Err(msg) => {
            let _ = writeln!(err, "ERROR: {msg}");
            1
        }
    }
}

fn run(o: &Opts) -> Result<Option<String>, String> {
    let bam = o.bam.as_deref().unwrap_or("the BAM");
    let from_bam = || -> Result<Option<String>, String> {
        match &o.bam {
            Some(path) => bam_basecaller(&bam_header(path).map_err(|e| e.0)?)
                .map_err(|e| format!("{bam}: {e}")),
            None => Ok(None),
        }
    };
    if o.basecaller_only {
        if o.bam.is_none() {
            return Err("--basecaller-only needs --bam".to_owned());
        }
        return from_bam();
    }
    let models_dir = o.models_dir.as_deref().unwrap_or("/opt/models");
    let available = || {
        let names = models_in(models_dir);
        if names.is_empty() {
            format!("no model directory found under {models_dir}")
        } else {
            format!("models in {models_dir}: {}", names.join(" "))
        }
    };
    let (path, why) = if let Some(model) = o.model.as_deref().filter(|m| !m.is_empty()) {
        (model.to_owned(), "set by CLAIR3_MODEL".to_owned())
    } else {
        let (basecaller, source) = if let Some(b) =
            o.basecaller.as_deref().filter(|b| !b.is_empty())
        {
            (b.to_owned(), "BASECALL_MODEL")
        } else {
            match from_bam() {
                Ok(Some(one)) => (one, "the BAM's read group"),
                Ok(None) => {
                    return Err(format!(
                        "{bam} names no basecaller model (no basecall_model= in an @RG DS field), so the Clair3 model cannot be chosen. Set BASECALL_MODEL to the Dorado model the reads were basecalled with (for example dna_r10.4.1_e8.2_400bps_sup@v5.0.0), or CLAIR3_MODEL to a model directory; {}.",
                        available()
                    ));
                }
                Err(e) => {
                    return Err(format!(
                        "{e}. Set CLAIR3_MODEL to the model directory to use, or BASECALL_MODEL; {}.",
                        available()
                    ));
                }
            }
        };
        let Some(dir) = clair3_model_for(&basecaller) else {
            return Err(format!(
                "Clair3 has no model for the basecaller model {basecaller} (from {source}). Clair3 models exist for R10.4.1 400 bps HAC and SUP and for R9.4.1 HAC and SUP basecalls, not for fast ones. Set CLAIR3_MODEL to a model directory if you know which fits; {}.",
                available()
            ));
        };
        (
            format!("{}/{dir}", models_dir.trim_end_matches('/')),
            format!("basecaller {basecaller}, from {source}"),
        )
    };
    if !Path::new(&path).is_dir() {
        return Err(format!(
            "the Clair3 model {path} ({why}) is not in the image. Set CLAIR3_MODEL to a model directory that exists; {}.",
            available()
        ));
    }
    Ok(Some(path))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(b: &str) -> Option<String> {
        clair3_model_for(b)
    }

    #[test]
    fn r10_sup_and_hac_map_by_version() {
        for (basecaller, want) in [
            (
                "dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
                "r1041_e82_400bps_sup_v410",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_hac@v4.1.0",
                "r1041_e82_400bps_hac_v410",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_sup@v4.3.0",
                "r1041_e82_400bps_sup_v430",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
                "r1041_e82_400bps_sup_v500",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
                "r1041_e82_400bps_hac_v500",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
                "r1041_e82_400bps_sup_v520",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
                "r1041_e82_400bps_hac_v520",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_hac@v6.0.0",
                "r1041_e82_400bps_hac_v600",
            ),
            (
                "dna_r10.4.1_e8.2_400bps_sup@v6.0.0",
                "r1041_e82_400bps_sup_v600",
            ),
            ("r10.4.1_e8.2_400bps_sup@v5.0", "r1041_e82_400bps_sup_v500"),
            (
                " dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n",
                "r1041_e82_400bps_sup_v500",
            ),
        ] {
            assert_eq!(model(basecaller).as_deref(), Some(want), "{basecaller}");
        }
    }

    #[test]
    fn r9_maps_to_the_two_r9_models() {
        for (basecaller, want) in [
            ("dna_r9.4.1_e8_sup@v3.3", "r941_prom_sup_g5014"),
            ("dna_r9.4.1_e8_sup@v3.6", "r941_prom_sup_g5014"),
            ("dna_r9.4.1_450bps_sup", "r941_prom_sup_g5014"),
            ("dna_r9.4.1_450bps_sup_prom", "r941_prom_sup_g5014"),
            ("dna_r9.4.1_e8_hac@v3.3", "r941_prom_hac_g360+g422"),
            ("dna_r9.4.1_450bps_hac_prom", "r941_prom_hac_g360+g422"),
        ] {
            assert_eq!(model(basecaller).as_deref(), Some(want), "{basecaller}");
        }
    }

    #[test]
    fn unknown_basecallers_map_to_nothing() {
        for basecaller in [
            "",
            "dna_r10.4.1_e8.2_400bps_fast@v5.0.0",
            "dna_r9.4.1_e8_fast@v3.4",
            "dna_r10.4.1_e8.2_260bps_sup@v4.0.0",
            "dna_r10.4.1_e8.2_400bps_sup",
            "dna_r10.4.1_e8.2_400bps_sup@v10.0.0",
            "dna_r10.4.1_e8.2_400bps_sup@5.0.0",
            "dna_r10.4.1_e8.2_400bps_sup@v5",
            "rna004_130bps_sup@v5.0.0",
            "dna_r10.4.1_e8.2_400bps_sup_v500",
            "hifi_revio",
        ] {
            assert_eq!(model(basecaller), None, "{basecaller:?}");
        }
    }

    #[test]
    fn basecallers_come_from_rg_ds_only() {
        let header = "@HD\tVN:1.6\n\
            @RG\tID:a\tPL:ONT\tDS:runid=x basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 modbase_models=m\n\
            @RG\tID:b\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0\r\n\
            @PG\tID:p\tDS:basecall_model=not_a_read_group\n\
            @CO\tbasecall_model=comment\n";
        assert_eq!(
            basecallers(header),
            (vec!["dna_r10.4.1_e8.2_400bps_sup@v5.0.0".to_owned()], 0)
        );
        assert_eq!(
            bam_basecaller(header),
            Ok(Some("dna_r10.4.1_e8.2_400bps_sup@v5.0.0".to_owned()))
        );
        let two = "@RG\tID:a\tDS:basecall_model=x@v1.0.0\n@RG\tID:b\tDS:basecall_model=y@v1.0.0\n";
        assert_eq!(basecallers(two).0.len(), 2);
        assert!(bam_basecaller(two).is_err());
        assert_eq!(basecallers("@RG\tID:a\tPL:ONT\n"), (Vec::new(), 1));
        assert_eq!(bam_basecaller("@RG\tID:a\tPL:ONT\n"), Ok(None));
        assert_eq!(bam_basecaller("@HD\tVN:1.6\n"), Ok(None));
        // One read group named, one not: the other's reads may be from any
        // basecaller, so no model is inferred.
        let partial = "@RG\tID:a\tDS:basecall_model=x@v1.0.0\n@RG\tID:b\tDS:runid=r2\n";
        assert_eq!(basecallers(partial), (vec!["x@v1.0.0".to_owned()], 1));
        let e = bam_basecaller(partial).unwrap_err();
        assert!(
            e.contains("1 read group(s) name no basecaller model"),
            "{e}"
        );
    }

    fn scratch(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("og-longread-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join("r1041_e82_400bps_sup_v500")).unwrap();
        std::fs::create_dir_all(d.join("hifi_revio")).unwrap();
        d
    }

    /// A BAM holding only a header: BGZF is gzip, which the reader takes.
    fn bam(dir: &Path, name: &str, header: &str) -> String {
        use std::io::Write as _;
        let mut raw = b"BAM\x01".to_vec();
        raw.extend_from_slice(&i32::try_from(header.len()).unwrap().to_le_bytes());
        raw.extend_from_slice(header.as_bytes());
        raw.extend_from_slice(&0i32.to_le_bytes());
        let path = dir.join(name);
        let mut gz = flate2::write::GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        gz.write_all(&raw).unwrap();
        gz.finish().unwrap();
        path.to_string_lossy().into_owned()
    }

    fn cli(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = clair3_model_cli(&args, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    #[test]
    fn the_cli_picks_checks_and_refuses() {
        let d = scratch("cli");
        let models = d.to_string_lossy().into_owned();
        let with = bam(
            &d,
            "with.bam",
            "@RG\tID:a\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n",
        );
        let (code, out, _) = cli(&["--bam", &with, "--models-dir", &models]);
        assert_eq!(code, 0);
        assert_eq!(out, format!("{models}/r1041_e82_400bps_sup_v500\n"));

        // BASECALL_MODEL wins over the BAM, and a missing directory stops.
        let (code, out, err) = cli(&[
            "--bam",
            &with,
            "--basecaller",
            "dna_r10.4.1_e8.2_400bps_hac@v4.1.0",
            "--models-dir",
            &models,
        ]);
        assert_eq!((code, out.as_str()), (1, ""));
        assert!(err.contains("r1041_e82_400bps_hac_v410"), "{err}");
        assert!(err.contains("is not in the image"), "{err}");
        assert!(
            err.contains("hifi_revio r1041_e82_400bps_sup_v500"),
            "{err}"
        );

        // CLAIR3_MODEL wins over both, and is checked too.
        let own = format!("{models}/hifi_revio");
        let (code, out, _) = cli(&[
            "--bam",
            &with,
            "--basecaller",
            "x",
            "--model",
            &own,
            "--models-dir",
            &models,
        ]);
        assert_eq!((code, out), (0, format!("{own}\n")));
        let (code, _, err) = cli(&["--model", "/no/such/model", "--models-dir", &models]);
        assert_eq!(code, 1);
        assert!(err.contains("CLAIR3_MODEL"), "{err}");

        // No basecaller anywhere, and a fast basecaller: both stop.
        let without = bam(&d, "without.bam", "@RG\tID:a\tPL:ONT\n");
        let (code, _, err) = cli(&["--bam", &without, "--models-dir", &models]);
        assert_eq!(code, 1);
        assert!(
            err.contains("BASECALL_MODEL") && err.contains("CLAIR3_MODEL"),
            "{err}"
        );
        let (code, _, err) = cli(&[
            "--basecaller",
            "dna_r10.4.1_e8.2_400bps_fast@v5.0.0",
            "--models-dir",
            &models,
        ]);
        assert_eq!(code, 1);
        assert!(err.contains("has no model"), "{err}");

        // --basecaller-only prints what the header names, or nothing.
        let (code, out, _) = cli(&["--bam", &with, "--basecaller-only"]);
        assert_eq!(
            (code, out.as_str()),
            (0, "dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n")
        );
        let (code, out, _) = cli(&["--bam", &without, "--basecaller-only"]);
        assert_eq!((code, out.as_str()), (0, ""));

        // A BAM with one read group named and one not: no model is inferred,
        // in either mode, but BASECALL_MODEL still decides.
        let partial = bam(
            &d,
            "partial.bam",
            "@RG\tID:a\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n@RG\tID:b\tPL:ONT\n",
        );
        let (code, out, err) = cli(&["--bam", &partial, "--models-dir", &models]);
        assert_eq!((code, out.as_str()), (1, ""));
        assert!(
            err.contains("name no basecaller model while others name"),
            "{err}"
        );
        let (code, out, _) = cli(&["--bam", &partial, "--basecaller-only"]);
        assert_eq!((code, out.as_str()), (1, ""));
        let (code, out, _) = cli(&[
            "--bam",
            &partial,
            "--basecaller",
            "dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
            "--models-dir",
            &models,
        ]);
        assert_eq!(
            (code, out),
            (0, format!("{models}/r1041_e82_400bps_sup_v500\n"))
        );

        // A usage error is exit 2.
        let (code, _, err) = cli(&["--og-scaffold-probe"]);
        assert_eq!(code, 2);
        assert!(err.contains("unrecognized arguments"), "{err}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
