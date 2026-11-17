//! `og facts`: the `ownhealth.genome-facts/v1` file of a finished run, the
//! one file OwnHealth imports (docs/ownhealth.md). OwnHealth owns the format;
//! `assets/ownhealth/genome-facts.v1.json` is its schema, copied unchanged.
//!
//! ```text
//! og facts --sample S --sample-dir DIR --reference FASTA [--assets DIR] [--out FILE]
//! ```
//!
//! Reads the sample folder's stable outputs ([`build`]) and writes
//! `DIR/facts/S_genome_facts.json` (or `--out`). The reference is the
//! uncompressed GRCh38 FASTA the run used, with its `.fai`: the VRS
//! identifiers are worked out on it. What is left out, and why, goes to
//! stderr. Exit 0 written, 1 a missing or unreadable input, 2 bad arguments.

pub mod build;
pub mod panels;
pub mod rules;
pub mod vrs;

use std::io::Write;
use std::path::{Path, PathBuf};

use og_intervals::fasta::Indexed;
use pycompat::json::{self, Dump};

const USAGE: &str =
    "usage: og facts --sample S --sample-dir DIR --reference FASTA [--assets DIR] [--out FILE]\n";

/// The uncompressed reference FASTA, read through its `.fai`.
pub struct Fasta(pub Indexed);

impl build::Genome for Fasta {
    fn get(&mut self, chrom: &str, beg: u64, end: u64) -> Result<Vec<u8>, String> {
        self.0
            .fetch(chrom, beg, end)
            .map_err(|e| e.to_string())?
            .ok_or_else(|| format!("{chrom} is not in the reference"))
    }
}

/// The reference must name GRCh38's chromosomes (`chr1` to `chrM`), and each
/// it has must have GRCh38's length, or every VRS identifier would name the
/// wrong sequence.
pub fn check_build(fasta: &Indexed, p: &panels::Panels) -> Result<(), String> {
    let mut found = 0;
    for (chrom, (_, len)) in &p.refget {
        if let Some(e) = fasta.entry(chrom) {
            if e.len != *len {
                return Err(format!(
                    "the reference's {chrom} has {} bases, GRCh38's has {len}: not a GRCh38 reference",
                    e.len
                ));
            }
            found += 1;
        }
    }
    if found == 0 {
        return Err(
            "the reference names none of GRCh38's chromosomes (chr1 to chr22, chrX, chrY, chrM): not a GRCh38 reference"
                .to_owned(),
        );
    }
    Ok(())
}

/// The `assets/` folder of the checkout the binary sits in, as `og summary
/// collect` finds it.
fn default_assets() -> PathBuf {
    let exe = std::env::current_exe().unwrap_or_default();
    exe.ancestors()
        .skip(1)
        .map(|d| d.join("assets"))
        .find(|p| p.join("ownhealth/genome-facts.v1.json").is_file())
        .unwrap_or_else(|| PathBuf::from("assets"))
}

struct Args {
    sample: String,
    dir: PathBuf,
    reference: PathBuf,
    assets: PathBuf,
    out: Option<PathBuf>,
}

fn parse(args: &[String]) -> Result<Args, String> {
    let (mut sample, mut dir, mut reference, mut assets, mut out) = (None, None, None, None, None);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let slot = match a.as_str() {
            "--sample" => &mut sample,
            "--sample-dir" => &mut dir,
            "--reference" => &mut reference,
            "--assets" => &mut assets,
            "--out" => &mut out,
            other => return Err(format!("unrecognized argument: {other}")),
        };
        *slot = Some(
            it.next()
                .ok_or_else(|| format!("{a} needs a value"))?
                .clone(),
        );
    }
    let need = |v: Option<String>, k: &str| v.ok_or_else(|| format!("{k} is required"));
    Ok(Args {
        sample: need(sample, "--sample")?,
        dir: need(dir, "--sample-dir")?.into(),
        reference: need(reference, "--reference")?.into(),
        assets: assets.map_or_else(default_assets, PathBuf::from),
        out: out.map(PathBuf::from),
    })
}

fn run(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> Result<PathBuf, String> {
    let p = panels::Panels::load(&a.assets)?;
    let fasta =
        Indexed::open(&a.reference).map_err(|e| format!("{}: {e}", a.reference.display()))?;
    check_build(&fasta, &p)?;
    let built = build::build(&a.sample, &a.dir, &p, &mut Fasta(fasta))?;
    for w in &built.warnings {
        let _ = writeln!(err, "og facts: {w}");
    }
    let target = a.out.clone().unwrap_or_else(|| {
        a.dir
            .join("facts")
            .join(format!("{}_genome_facts.json", a.sample))
    });
    let body = json::dumps(
        &built.doc,
        &Dump {
            indent: Some(1),
            ensure_ascii: false,
            ..Dump::default()
        },
    )
    .map_err(|e| format!("{e:?}"))?
        + "\n";
    write_atomic(&target, body.as_bytes())?;
    let n = match &built.doc {
        json::Value::Object(m) => m.get("facts").map_or(0, |f| match f {
            json::Value::Array(a) => a.len(),
            _ => 0,
        }),
        _ => 0,
    };
    let _ = writeln!(
        out,
        "Genome facts for OwnHealth written: {} ({n} facts)",
        target.display()
    );
    Ok(target)
}

/// Writes `path` through a temporary file of this process beside it, so a
/// reader never sees half a file and two writers never share one.
fn write_atomic(path: &Path, body: &[u8]) -> Result<(), String> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d).map_err(|e| format!("{}: {e}", d.display()))?;
    }
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(format!(".{}.part", std::process::id()));
    let tmp = PathBuf::from(tmp);
    let written = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&tmp)
        .and_then(|mut f| f.write_all(body))
        .and_then(|()| std::fs::rename(&tmp, path));
    written.map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        format!("{}: {e}", path.display())
    })
}

/// Runs `og facts ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match parse(args) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og facts: error: {e}");
            return 2;
        }
    };
    match run(&a, out, err) {
        Ok(_) => 0,
        Err(e) => {
            let _ = writeln!(err, "og facts: ERROR: {e}");
            1
        }
    }
}
