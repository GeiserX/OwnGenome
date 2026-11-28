//! `og setup install [FLAG] <genome_dir>`: setup.sh's installers, run by og.
//!
//! It takes exactly setup.sh's arguments and parses them the same way
//! ([`parse`]), with one difference: a second flag after the first is
//! refused instead of being read as the install folder. A flag og has
//! ported ([`ported`]) runs here; any other is handed to
//! `scripts/setup.sh` with `OG_SETUP=bash`, so `og setup install` takes
//! every flag. setup.sh hands a ported flag to og the same way
//! (`OG_SETUP=bash` keeps it in bash: tests/setup/ab.sh compares the two).
//!
//! Before an installer runs, og sets up what scripts/lib/common.sh sets up
//! for setup.sh ([`Run::new`]): umask 077, the install folder without its
//! trailing slash, the downloads folder (OG_DOWNLOADS, made absolute), the
//! reference path and its container path, and setup.sh's check that the
//! kept downloads fit (`og setup plan --check`).
//!
//! Exit status: an installer's own (0 installed, 1 failed); 1 for the usage
//! text (no install folder) or a bad `--refresh`; 2 for an unknown flag, a
//! second flag, a bad SETUP_KEEP_DOWNLOADS or OG_SETUP, or a reference
//! outside the install folder; 3 is a full setup's partial status
//! ([`Missing`]).
//!
//! Each later slice of the port fills one module, declared here once:
//! [`clinical`], [`annotation`], [`panels`], [`tools`] and [`base`] (the
//! default run, `--everything`, `--refresh clinvar`). [`small`] holds the
//! eight small installers. A module lists the flags it runs in its `PORTED`
//! and answers them in its `run`, so adding one edits no shared list.

pub mod annotation;
pub mod base;
pub mod clinical;
pub mod panels;
pub mod small;
pub mod tools;

use crate::components;
use crate::fetch::{
    self, Ctx, Env as FetchEnv, Fetch, Kind, Knobs, OnWrong, Sources, in_downloads, source_for,
    versions_pins,
};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The flags of setup.sh's opt-in case, in its order.
pub const OPT_IN: [&str; 25] = [
    "cyrius",
    "prs-academic",
    "parascopy-data",
    "kir-data",
    "ancestry-panel",
    "ancestry-fine-panel",
    "yleaf-data",
    "consumer-sites",
    "ancient-panel",
    "viral-data",
    "masked-reference",
    "exomiser-data",
    "chip-data",
    "strchive-data",
    "megane-kmers",
    "mocha-data",
    "locityper-data",
    "vep-cache",
    "annotsv-data",
    "annotation-scores",
    "somatic-data",
    "bwa-index",
    "cnvpytor-data",
    "pypgx-data",
    "archaic",
];

/// What one setup run does, from its first argument.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mode {
    /// `setup.sh <dir>`.
    Full,
    Everything,
    NoClinical,
    /// One opt-in installer (`--clinical-data` included).
    OptIn(String),
    PullOnly,
    SampleQc,
    /// `--refresh clinvar`.
    Refresh,
    /// `--og`.
    Og,
}

impl Mode {
    /// The component key: a flag without its dashes, `setup` for a plain run.
    pub fn key(&self) -> &str {
        match self {
            Mode::Full => "setup",
            Mode::Everything => "everything",
            Mode::NoClinical => "no-clinical-data",
            Mode::OptIn(k) => k,
            Mode::PullOnly => "pull-only",
            Mode::SampleQc => "sample-qc-data",
            Mode::Refresh => "refresh",
            Mode::Og => "og",
        }
    }

    /// The flag that starts it ("" for a plain run).
    pub fn flag(&self) -> String {
        match self {
            Mode::Full => String::new(),
            Mode::Refresh => "--refresh".to_owned(),
            m => format!("--{}", m.key()),
        }
    }
}

/// setup.sh's arguments, parsed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Parsed {
    pub mode: Mode,
    /// The install folder as given ("" for `--pull-only` without one).
    pub genome_dir: String,
    /// SETUP_KEEP_DOWNLOADS after `--no-keep-downloads`.
    pub keep_downloads: bool,
    /// The arguments setup.sh hands on: the flag, then the rest.
    pub argv: Vec<String>,
}

/// Why parsing stopped: the exit code, and what to print.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Stop {
    pub code: u8,
    pub stdout: String,
    pub stderr: String,
}

/// setup.sh's usage text, as its `echo` lines hold it: `$0` is the script,
/// `${NAME:-}` and `${NAME}` are read from versions.env.
const USAGE: &[&str] = &[
    "Usage: $0 <genome_dir>",
    "       $0 --everything <genome_dir>",
    "       $0 --no-clinical-data <genome_dir>",
    "       $0 --clinical-data <genome_dir>",
    "       $0 --pull-only",
    "       $0 --og",
    "       $0 --refresh clinvar <genome_dir>",
    "       $0 --sample-qc-data <genome_dir>",
    "       $0 --cyrius | --prs-academic | --parascopy-data | --kir-data | --ancestry-panel | --ancestry-fine-panel | --yleaf-data | --consumer-sites | --ancient-panel | --viral-data | --masked-reference | --exomiser-data | --chip-data | --strchive-data | --megane-kmers | --mocha-data | --locityper-data | --archaic | --vep-cache | --annotsv-data <genome_dir>",
    "",
    "  <genome_dir>  Where to store reference data and sample outputs.",
    "                Needs at least 500 GB free space per sample.",
    "  --everything  Run full setup and every installer in assets/everything.tsv.",
    "  --annotation-scores  Install CADD, masked SpliceAI, REVEL and AlphaMissense (~175 GiB).",
    "  --somatic-data       Install GATK somatic resources for steps 29 and 49.",
    "  --bwa-index          Build the classic BWA index for GRIDSS.",
    "  --cnvpytor-data      Install CNVpytor reference files (step 18).",
    "  --pypgx-data         Install the upstream PyPGx bundle (step 32).",
    "  --no-clinical-data  Skip both VEP caches and the PCGR/CPSR bundle.",
    "  --no-keep-downloads  Delete each raw download once it is installed (SETUP_KEEP_DOWNLOADS=false).",
    "                Without it every raw download stays in the downloads folder, OG_DOWNLOADS",
    "                (default <genome_dir>/downloads), so a rerun downloads nothing again; this keeps",
    "                about 300 GB more for --everything. bin/og setup plan lists what a setup downloads,",
    "                and bin/og setup prune deletes the kept downloads. Goes with any other flag.",
    "  --clinical-data     Install only those three downloads for an existing setup.",
    "  --pull-only   Pull the Docker images listed in versions.env and exit.",
    "  --og          Put OwnGenome's og binaries in bin/.og, for Linux and for a macOS live-site host",
    "                (built from rust/, or a release's) and exit.",
    "  --refresh clinvar",
    "                Download NCBI's current ClinVar, rebuild the files made from it",
    "                and record its release date in <genome_dir>/clinvar/RELEASE.",
    "  --sample-qc-data",
    "                Install only somalier's sites and VerifyBamID2's panel (step 33) and exit.",
    "  --cyrius      Install Cyrius for the opt-in step 21 (PyPI, hash-locked) and exit.",
    "                Cyrius is under the PolyForm Strict licence: non-commercial use only.",
    "  --prs-academic",
    "                Install the four PGS Catalog scores of step 25 that are free for academic research use only",
    "                (PGS000014, PGS000016, PGS000017, PGS000027; PRS_ACADEMIC=true scores them) and exit.",
    "  --parascopy-data",
    "                Install Parascopy's GRCh38 homology table and models (step 35, ~50 MB), pull its image and exit.",
    "  --kir-data    Install the IPD-KIR ${KIR_DB_RELEASE:-} database (KIR=true in step 08, ~40 MB) and exit.",
    "  --ancestry-panel",
    "                Install pgsc_calc's ancestry reference panel ${PGSC_PANEL:-} (~7 GB download): step 26, and",
    "                percentiles instead of raw scores in step 25; then exit.",
    "  --ancestry-fine-panel",
    "                Install the reference of step 26b, gnomAD's HGDP + 1000 Genomes panel (~19 GB is read, a few GB",
    "                stored, hours on a home connection), and its images; then exit.",
    "  --yleaf-data  Install Yleaf ${YLEAF_DATA_VERSION:-}'s Y marker tables and trees and the current YFull tree (step 37, about 210 MB to download, 220 MB on disk), pull its image and exit.",
    "  --consumer-sites",
    "                Build the site table of the raw-data export (step 40): the 23andMe v5 site list joined with",
    "                dbSNP ${CONSUMER_DBSNP_BUILD:-}'s GRCh38 VCF (~30 GB download, kept in the downloads folder), with UCSC's",
    "                hg19ToHg38 chain for each site's strand; then exit.",
    "  --archaic     Build the archaic table of step 43: the Altai and Vindija Neanderthals and the Altai Denisovan",
    "                (MPI-EVA, ~165 GB, kept in the downloads folder) with 1000 Genomes frequencies, lifted to GRCh38",
    "                (hours on a home connection; --no-keep-downloads keeps 13 GB at most at a time); then exit.",
    "  --ancient-panel",
    "                Build the panel of step 44 from the Allen Ancient DNA Resource ${AADR_VERSION:-} (1240K, ~7.3 GB download,",
    "                kept in the downloads folder; needs about 12 GB of free RAM) and pull EIGENSOFT; then exit.",
    "  --viral-data  Install the human viral genomes and phiX of step 47 (NCBI RefSeq, VIRAL_REF_ACCESSIONS, ~2 MB) and exit.",
    "  --masked-reference",
    "                Build the default reference with the false duplications of GRCh38 masked (KCNE1, CBS, U2AF1, CRYAA",
    "                and others become readable; ~3.2 GB beside it, a few minutes) and exit. Samples must be realigned to it.",
    "  --exomiser-data",
    "                Install Exomiser ${EXOMISER_DATA_VERSION:-}'s hg38 and phenotype data (step 48: about 37 GB to download, 58 GB",
    "                unpacked, about 72 GB free while it unpacks), pull its image and exit.",
    "  --chip-data   Install Ensembl ${ENSEMBL_GFF_RELEASE:-}'s GRCh38 gene annotation (step 49, the CHIP genes' MANE Select",
    "                transcripts; ~110 MB) and exit.",
    "  --strchive-data",
    "                Install only the STRchive ${STRCHIVE_VERSION:-} repeat catalog of steps 09 and 9b (~0.5 MB) and exit.",
    "  --megane-kmers",
    "                Build MEGAnE's k-mer set of the reference for step 50 (Alu, L1 and SVA insertions; ~1.8 GB,",
    "                a few minutes and about 5 GB of RAM, with og), pull MEGAnE's image (~6.4 GB) and exit.",
    "  --mocha-data  Install MoChA's phasing panel for step 51 (1000 Genomes high coverage, ${MOCHA_PANEL_VERSION:-}: ~30 GB read,",
    "                one chromosome at a time) and the genetic map, pull its images and exit.",
    "  --locityper-data",
    "                Build Locityper's database for step 53: the k-mer counts of the reference (several GB of RAM",
    "                once, a few GB on disk), the HPRC ${LOCITYPER_PANEL_VERSION:-} haplotypes of the loci in",
    "                assets/locityper_loci.bed (read over https, tens of MB), their pairwise alignments; pull its image and exit.",
    "  --annotsv-data",
    "                Install AnnotSV's data for step 5 and exit: the bundle when it is missing (~5.3 GB), then",
    "                gnomAD-SV ${GNOMAD_SV_VERSION:-}'s common SVs (~1.7 GB download) and MANE ${MANE_RELEASE:-}'s Select transcripts.",
    "  --vep-cache   Install the VEP ${VEP_CACHE_RELEASE:-} cache (step 13 and the Nextflow VEP: ~26 GB download, ~30 GB",
    "                unpacked, both on disk during the install) and exit. Without it a run skips VEP, vcfanno,",
    "                the clinical filter and slivar.",
    "",
    "Download and use terms (see docs/licences.md for the sources):",
    "  ${ANNOTSV_TERMS}",
    "  PCGR/CPSR reference data includes third-party databases with their own terms.",
    "  ${EXOMISER_TERMS}",
    "  The REVEL authors make the scores freely available for non-commercial use.",
    "  ${UCSC_CHAIN_TERMS}",
    "  Illumina's DRAGEN ORA Decompression Software requires acceptance of its end-user licence.",
    "  The PacBio Software License Agreement permits Paraphase use only with data generated on a PacBio instrument or provided by PacBio.",
    "",
    "Example:",
    "  ./scripts/setup.sh /data/genomics",
    "  ./scripts/setup.sh ~/genome_data",
];

/// The usage text of SCRIPT, with versions.env's values from VERSIONS.
/// Only the template's own `$0` and `${NAME}` are expanded: a value that
/// holds `${...}` is written as it is.
pub fn usage(script: &str, versions: &dyn Fn(&str) -> String) -> String {
    let mut s = String::new();
    for line in USAGE {
        let mut rest = *line;
        while let Some(i) = rest.find('$') {
            s.push_str(&rest[..i]);
            let at = &rest[i..];
            if let Some(after) = at.strip_prefix("$0") {
                s.push_str(script);
                rest = after;
            } else if let Some(j) = at.strip_prefix("${").and_then(|a| a.find('}')) {
                let inner = &at[2..2 + j];
                let name = inner.split(":-").next().unwrap_or(inner);
                s.push_str(&match name {
                    "ANNOTSV_TERMS" => components::ANNOTSV_TERMS.to_owned(),
                    "EXOMISER_TERMS" => components::EXOMISER_TERMS.to_owned(),
                    "UCSC_CHAIN_TERMS" => components::UCSC_CHAIN_TERMS.to_owned(),
                    n => versions(n),
                });
                rest = &at[3 + j..];
            } else {
                s.push('$');
                rest = &at[1..];
            }
        }
        s.push_str(rest);
        s.push('\n');
    }
    s
}

/// Parses setup.sh's arguments as setup.sh does. ENV reads
/// SETUP_KEEP_DOWNLOADS and GENOME_DIR; VERSIONS the values the usage text
/// names; SCRIPT is setup.sh's path as it was called (`$0`).
pub fn parse(
    args: &[String],
    env: &dyn Fn(&str) -> Option<String>,
    versions: &dyn Fn(&str) -> String,
    script: &str,
) -> Result<Parsed, Stop> {
    let stop = |code: u8, stdout: String, stderr: String| Stop {
        code,
        stdout,
        stderr,
    };
    let mut keep = env("SETUP_KEEP_DOWNLOADS")
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| "true".to_owned());
    let mut rest: Vec<String> = Vec::new();
    for a in args {
        if a == "--no-keep-downloads" {
            keep = "false".to_owned();
        } else {
            rest.push(a.clone());
        }
    }
    if keep != "true" && keep != "false" {
        return Err(stop(
            2,
            String::new(),
            format!("ERROR: SETUP_KEEP_DOWNLOADS must be true or false, got '{keep}'\n"),
        ));
    }
    let first = rest.first().cloned().unwrap_or_default();
    let mut shift = 1;
    let mode = match first.as_str() {
        "--everything" => Mode::Everything,
        "--og" => Mode::Og,
        f if f.strip_prefix("--").is_some_and(|k| OPT_IN.contains(&k)) => {
            Mode::OptIn(f[2..].to_owned())
        }
        "--clinical-data" => Mode::OptIn("clinical-data".to_owned()),
        "--no-clinical-data" => Mode::NoClinical,
        "--pull-only" => Mode::PullOnly,
        "--sample-qc-data" => Mode::SampleQc,
        "--refresh" => {
            if rest.get(1).map(String::as_str) != Some("clinvar") {
                return Err(stop(
                    1,
                    String::new(),
                    format!(
                        "Usage: {script} --refresh clinvar <genome_dir>\n  (clinvar is the one database this flag refreshes)\n"
                    ),
                ));
            }
            shift = 2;
            Mode::Refresh
        }
        f if f.starts_with("--") => {
            return Err(stop(
                2,
                String::new(),
                format!("ERROR: unknown setup flag: {f}\n"),
            ));
        }
        _ => {
            shift = 0;
            Mode::Full
        }
    };
    if mode == Mode::Og {
        return Ok(Parsed {
            mode,
            genome_dir: String::new(),
            keep_downloads: keep == "true",
            argv: rest,
        });
    }
    let after: Vec<String> = rest.iter().skip(shift).cloned().collect();
    // One flag at a time: setup.sh would read a second one as the folder.
    // A flag og hands back to setup.sh keeps setup.sh's reading.
    if ported(mode.key())
        && let Some(second) = after.iter().find(|a| a.starts_with("--"))
    {
        let by = if first.starts_with("--") {
            first.clone()
        } else {
            format!("the install folder {first}")
        };
        return Err(stop(
            2,
            String::new(),
            format!(
                "ERROR: setup takes one flag at a time: {second} comes after {by}. Run setup once for each.\n"
            ),
        ));
    }
    let genome_dir = after
        .first()
        .cloned()
        .filter(|g| !g.is_empty())
        .or_else(|| env("GENOME_DIR").filter(|g| !g.is_empty()))
        .unwrap_or_default();
    if genome_dir.is_empty() && mode != Mode::PullOnly {
        return Err(stop(1, usage(script, versions), String::new()));
    }
    Ok(Parsed {
        mode,
        genome_dir,
        keep_downloads: keep == "true",
        argv: rest,
    })
}

/// Whether og runs KEY (a component key) itself.
pub fn ported(key: &str) -> bool {
    [
        small::PORTED,
        clinical::PORTED,
        annotation::PORTED,
        panels::PORTED,
        tools::PORTED,
        base::PORTED,
    ]
    .iter()
    .any(|p| p.contains(&key))
}

/// A word as bash's `printf %q` writes it (bash 5.2, in a UTF-8 locale).
pub fn shell_quote(s: &str) -> String {
    if s.is_empty() {
        return "''".to_owned();
    }
    if s.chars().any(|c| c.is_control()) {
        let mut o = String::from("$'");
        for c in s.chars() {
            match c {
                '\x07' => o.push_str("\\a"),
                '\x08' => o.push_str("\\b"),
                '\t' => o.push_str("\\t"),
                '\n' => o.push_str("\\n"),
                '\x0b' => o.push_str("\\v"),
                '\x0c' => o.push_str("\\f"),
                '\r' => o.push_str("\\r"),
                '\x1b' => o.push_str("\\E"),
                '\'' => o.push_str("\\'"),
                '\\' => o.push_str("\\\\"),
                // Each byte of its UTF-8 form: U+0085 is $'\302\205',
                // which $'...' turns back into those two bytes.
                c if c.is_control() => {
                    for b in c.to_string().bytes() {
                        o.push_str(&format!("\\{b:03o}"));
                    }
                }
                c => o.push(c),
            }
        }
        o.push('\'');
        return o;
    }
    let mut o = String::new();
    for (i, c) in s.chars().enumerate() {
        let special = match c {
            'a'..='z' | 'A'..='Z' | '0'..='9' => false,
            '_' | '.' | '/' | ',' | ':' | '@' | '%' | '+' | '=' | '-' => false,
            '~' | '#' => i == 0,
            c if !c.is_ascii() => false,
            _ => true,
        };
        if special {
            o.push('\\');
        }
        o.push(c);
    }
    o
}

/// The environment of one install: the process's, with what common.sh
/// would have set (GENOME_DIR, OG_DOWNLOADS, SETUP_KEEP_DOWNLOADS, REF_FASTA).
#[derive(Debug, Clone, Default)]
pub struct Vars {
    map: BTreeMap<String, String>,
}

impl Vars {
    /// The process's environment. A variable whose name or value is not
    /// UTF-8 is left out (none of setup's settings can be one), where
    /// `std::env::vars` would panic.
    pub fn from_process() -> Vars {
        Vars {
            map: std::env::vars_os()
                .filter_map(|(k, v)| Some((k.into_string().ok()?, v.into_string().ok()?)))
                .collect(),
        }
    }

    pub fn from_map(map: BTreeMap<String, String>) -> Vars {
        Vars { map }
    }

    /// The value, set or empty.
    pub fn get(&self, k: &str) -> Option<&str> {
        self.map.get(k).map(String::as_str)
    }

    /// The value when it is not empty (`${NAME:-}`).
    pub fn nonempty(&self, k: &str) -> Option<&str> {
        self.get(k).filter(|v| !v.is_empty())
    }

    pub fn set(&mut self, k: &str, v: &str) {
        self.map.insert(k.to_owned(), v.to_owned());
    }

    /// A lookup for the fetch and plan code.
    pub fn lookup(&self) -> impl Fn(&str) -> Option<String> + 'static {
        let m = self.map.clone();
        move |k| m.get(k).cloned()
    }

    pub fn boxed(&self) -> FetchEnv {
        Box::new(self.lookup())
    }
}

/// Everything an installer reads: the checkout, the environment and
/// versions.env.
pub struct Run {
    pub root: PathBuf,
    pub vars: Vars,
    /// setup.sh as it was called, for messages that name it (`$0`).
    pub script: String,
    /// The install folder, without its trailing slash ("" for none).
    pub genome: String,
    pub pins: indexmap::IndexMap<String, String>,
    /// The og binary, for the transfers setup ran as `og setup fetch --get`.
    pub og: PathBuf,
}

impl Run {
    /// What common.sh does once sourced, and setup.sh before its
    /// installers: umask 077, GENOME_DIR without a trailing slash, the
    /// downloads folder, REF_FASTA and its container path. `Err` carries the
    /// exit code after the message is written.
    pub fn new(
        root: &Path,
        mut vars: Vars,
        script: &str,
        genome_dir: &str,
        keep_downloads: bool,
        err: &mut dyn Write,
    ) -> Result<Run, u8> {
        rustix::process::umask(rustix::fs::Mode::from_raw_mode(0o077));
        let versions = root.join("versions.env");
        let pins = std::fs::read_to_string(&versions)
            .map(|t| versions_pins(&t))
            .unwrap_or_default();
        vars.set(
            "SETUP_KEEP_DOWNLOADS",
            if keep_downloads { "true" } else { "false" },
        );
        vars.set("PGP_ROOT", &root.to_string_lossy());
        let mut genome = String::new();
        if !genome_dir.is_empty() {
            genome = genome_dir
                .strip_suffix('/')
                .unwrap_or(genome_dir)
                .to_owned();
            vars.set("GENOME_DIR", &genome);
            let pwd = vars
                .nonempty("PWD")
                .map(str::to_owned)
                .or_else(|| {
                    std::env::current_dir()
                        .ok()
                        .map(|p| p.to_string_lossy().into_owned())
                })
                .unwrap_or_default();
            let dl = match vars.nonempty("OG_DOWNLOADS") {
                Some(d) if d.starts_with('/') => d.to_owned(),
                Some(d) => format!("{pwd}/{d}"),
                None => format!("{genome}/downloads"),
            };
            let dl = dl.strip_suffix('/').unwrap_or(&dl).to_owned();
            vars.set("OG_DOWNLOADS", &dl);
            let reference = match vars.nonempty("REF_FASTA") {
                Some(r) if r.starts_with('/') => r.to_owned(),
                Some(r) => format!("{genome}/{r}"),
                None => format!("{genome}/reference/GRCh38_no_alt_analysis_set.fasta"),
            };
            vars.set("REF_FASTA", &reference);
            if let Err(m) = crate::engine::cpath(&genome, Some(&dl), &reference) {
                let _ = writeln!(err, "{m}");
                return Err(2);
            }
        }
        let og = std::env::current_exe().unwrap_or_else(|_| root.join("bin/og"));
        Ok(Run {
            root: root.to_path_buf(),
            vars,
            script: script.to_owned(),
            genome,
            pins,
            og,
        })
    }

    /// A versions.env value, else the environment's (`${NAME:-}` after
    /// sourcing versions.env), else "".
    pub fn version(&self, k: &str) -> String {
        self.pins
            .get(k)
            .cloned()
            .or_else(|| self.vars.get(k).map(str::to_owned))
            .unwrap_or_default()
    }

    pub fn ctx(&self) -> Result<Ctx, String> {
        Ctx::new(&self.root, self.vars.boxed())
    }

    /// `source_of ID FIELD` (path or url) of common.sh; its error goes to ERR.
    pub fn source_of(&self, id: &str, field: &str, err: &mut dyn Write) -> Option<String> {
        let r = (|| -> Result<String, String> {
            let ctx = self.ctx()?;
            let s = Sources::read(&ctx.root)?;
            let row = s.get(id)?;
            match field {
                "path" => Ok(ctx.path(row, None)?.to_string_lossy().into_owned()),
                "url" => ctx.url(row, None),
                f => Err(format!("ERROR: unknown field {f}")),
            }
        })();
        match r {
            Ok(v) => Some(v),
            Err(m) => {
                let _ = writeln!(err, "{m}");
                None
            }
        }
    }

    /// `fetch URL DEST KIND VALUE` of common.sh (`og setup fetch`), in
    /// process: the same request, the same messages on ERR.
    pub fn fetch(
        &self,
        url: &str,
        dest: &Path,
        kind: &str,
        want: &str,
        err: &mut dyn Write,
    ) -> bool {
        let env = self.vars.lookup();
        let knobs = match Knobs::from(&env) {
            Ok(k) => k,
            Err(m) => {
                let _ = writeln!(err, "{m}");
                return false;
            }
        };
        let kind = match Kind::parse(kind) {
            Ok(k) => k,
            Err(m) => {
                let _ = writeln!(err, "{m}");
                return false;
            }
        };
        let on_wrong = if in_downloads(dest, &env) {
            OnWrong::Refuse
        } else {
            OnWrong::Replace
        };
        let row = self
            .ctx()
            .ok()
            .and_then(|c| Sources::read(&c.root).ok().map(|s| (c, s)))
            .and_then(|(c, s)| source_for(&c, &s, dest));
        let r = Fetch {
            url,
            dest,
            kind,
            want,
            on_wrong,
            source: row.as_ref().map(|(l, _)| l.as_str()),
            prune: row.as_ref().and_then(|(_, p)| p.as_deref()),
        };
        fetch::fetch(&r, &knobs, err).is_ok()
    }

    /// `_get URL -` of common.sh with its errors discarded: the body, or
    /// `None` when the transfer failed.
    pub fn get_quiet(&self, url: &str) -> Option<Vec<u8>> {
        let o = crate::child::output(
            Command::new(&self.og)
                .args(["setup", "fetch", "--get", url, "-"])
                .env("PGP_ROOT", &self.root)
                .stdin(std::process::Stdio::null())
                .stderr(std::process::Stdio::null()),
        )
        .ok()?;
        o.status.success().then_some(o.stdout)
    }
}

/// `[ -s PATH ]`: there, and not empty.
pub fn nonempty(p: &Path) -> bool {
    std::fs::metadata(p).is_ok_and(|m| m.len() > 0)
}

/// `_digest sha256 FILE`, "" when it cannot be read.
pub fn sha256(p: &Path) -> String {
    fetch::digest(Kind::Sha256, p).unwrap_or_default()
}

/// What a full setup could not install, one entry each with the command
/// that retries it (setup.sh's MISSING and print_missing).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Missing {
    items: Vec<(String, String)>,
}

/// The exit status of a full setup that left items missing.
pub const SETUP_PARTIAL_EXIT: u8 = 3;

impl Missing {
    pub fn note(&mut self, what: &str, retry: &str) {
        self.items.push((what.to_owned(), retry.to_owned()));
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// print_missing's text.
    pub fn text(&self) -> String {
        let mut s = String::from(
            "\nMissing (the steps named run without it or say it is missing; the warnings above say why):\n",
        );
        for (what, retry) in &self.items {
            s.push_str(&format!("  - {what}\n      retry: {retry}\n"));
        }
        s
    }
}

const HELP: &str = "usage: og setup install [FLAG] <genome_dir>

  setup.sh's installers, run by og: the same arguments as scripts/setup.sh
  (one flag, then the install folder; --no-keep-downloads anywhere). A
  flag og has ported runs here (og setup ported FLAG says which; og setup
  flags lists them all); any other is run by scripts/setup.sh with
  OG_SETUP=bash. OG_SETUP=bash runs every flag with setup.sh's bash.
";

/// `og setup install ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if matches!(args, [h] if h == "--help" || h == "-h") {
        let _ = write!(out, "{HELP}");
        return 0;
    }
    let root = fetch::root_of(None);
    let vars = Vars::from_process();
    let script = vars
        .nonempty("OG_SETUP_SCRIPT")
        .map(str::to_owned)
        .unwrap_or_else(|| root.join("scripts/setup.sh").to_string_lossy().into_owned());
    let bash = match vars.get("OG_SETUP").unwrap_or("") {
        "" | "og" => false,
        "bash" => true,
        other => {
            let _ = writeln!(
                err,
                "ERROR: OG_SETUP must be og or bash, got '{other}' (bash runs setup.sh's own installers)"
            );
            return 2;
        }
    };
    if bash {
        return hand_to_bash(&root, args, err);
    }
    let pins = std::fs::read_to_string(root.join("versions.env"))
        .map(|t| versions_pins(&t))
        .unwrap_or_default();
    let versions = |k: &str| {
        pins.get(k)
            .cloned()
            .or_else(|| std::env::var(k).ok())
            .unwrap_or_default()
    };
    let env = vars.lookup();
    let p = match parse(args, &env, &versions, &script) {
        Ok(p) => p,
        Err(s) => {
            let _ = write!(out, "{}", s.stdout);
            let _ = write!(err, "{}", s.stderr);
            return s.code;
        }
    };
    if !ported(p.mode.key()) {
        return hand_to_bash(&root, args, err);
    }
    let run = match Run::new(&root, vars, &script, &p.genome_dir, p.keep_downloads, err) {
        Ok(r) => r,
        Err(code) => return code,
    };
    if p.mode != Mode::PullOnly {
        // setup.sh builds og first, or takes OG; og is running, so only a
        // bad OG is left to refuse.
        if let Some(o) = run.vars.nonempty("OG")
            && !is_executable(Path::new(o))
        {
            let _ = writeln!(err, "ERROR: OG={o} is not an executable og");
            return 1;
        }
        if let Some(code) = disk_check(&run, &p, err) {
            return code;
        }
    }
    let key = p.mode.key().to_owned();
    [
        small::run,
        clinical::run,
        annotation::run,
        panels::run,
        tools::run,
        base::run,
    ]
    .iter()
    .find_map(|f| f(&key, &run, out, err))
    .unwrap_or_else(|| {
        let _ = writeln!(
            err,
            "ERROR: og lists {key} as ported but runs no installer for it"
        );
        2
    })
}

/// setup.sh's check that the kept downloads fit (`og setup plan --check`),
/// under the same conditions; `Some` exit code when it refuses.
fn disk_check(run: &Run, p: &Parsed, err: &mut dyn Write) -> Option<u8> {
    let v = &run.vars;
    if !p.keep_downloads
        || v.get("SETUP_DISK_CHECK") == Some("false")
        || p.mode == Mode::Refresh
        || v.get("SETUP_EVERYTHING_SCOPE") == Some("row")
    {
        return None;
    }
    let mut args: Vec<String> = match &p.mode {
        Mode::Everything => vec!["--all".into()],
        Mode::OptIn(k) => vec!["--for".into(), k.clone()],
        Mode::SampleQc => vec!["--for".into(), "sample-qc-data".into()],
        Mode::Full => vec!["--for".into(), "setup".into()],
        Mode::NoClinical => vec![
            "--for".into(),
            "setup".into(),
            "--skip".into(),
            "clinical-data".into(),
        ],
        _ => return None,
    };
    args.push("--check".into());
    args.extend(["--root".into(), run.root.to_string_lossy().into_owned()]);
    let code = crate::plan::cli_with(&args, v.boxed(), &mut std::io::sink(), err);
    (code != 0).then_some(1)
}

/// Runs scripts/setup.sh with OG_SETUP=bash and ARGS; its exit status.
/// OG names this og unless it is set already, so setup.sh uses the og that
/// runs it (a dev build too) instead of checking or rebuilding bin/og.
fn hand_to_bash(root: &Path, args: &[String], err: &mut dyn Write) -> u8 {
    let _ = std::io::stdout().flush();
    let script = root.join("scripts/setup.sh");
    let mut c = Command::new(&script);
    c.args(args).env("OG_SETUP", "bash");
    if std::env::var_os("OG").is_none_or(|o| o.is_empty())
        && let Ok(me) = std::env::current_exe()
    {
        c.env("OG", me);
    }
    match crate::child::status(&mut c) {
        Ok(s) => s.code().map_or(1, |c| u8::try_from(c).unwrap_or(1)),
        Err(e) => {
            let _ = writeln!(err, "ERROR: cannot run {}: {e}", script.display());
            1
        }
    }
}

/// `[ -x PATH ]` for a file.
pub fn is_executable(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    std::fs::metadata(p).is_ok_and(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
}

/// `command -v NAME`: a path is taken as it is, a name is looked up on PATH.
pub fn command_v(name: &str) -> bool {
    if name.contains('/') {
        is_executable(Path::new(name))
    } else {
        fetch::on_path(name).is_some()
    }
}

/// `og setup ported FLAG`: exit 0 when og runs FLAG itself, 1 when
/// setup.sh's bash does, 2 when FLAG is not a flag of setup.
pub fn ported_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let [flag] = args else {
        let _ = writeln!(err, "usage: og setup ported FLAG");
        return 2;
    };
    match components::by_flag(flag) {
        Some(c) if c.ported() => {
            let _ = writeln!(out, "{flag}: og");
            0
        }
        Some(_) => {
            let _ = writeln!(out, "{flag}: bash (scripts/setup.sh)");
            1
        }
        None => {
            let _ = writeln!(err, "ERROR: {flag} is not a flag of setup");
            2
        }
    }
}

/// `og setup flags [--long]`: every flag of setup, one per line; with
/// --long, a tab-separated row each: flag, who runs it (og or bash), steps,
/// images, sources, and the installed paths (relative to the install folder).
pub fn flags_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let long = match args {
        [] => false,
        [l] if l == "--long" => true,
        _ => {
            let _ = writeln!(err, "usage: og setup flags [--long]");
            return 2;
        }
    };
    if !long {
        for f in components::flags() {
            let _ = writeln!(out, "{f}");
        }
        return 0;
    }
    let root = fetch::root_of(None);
    let s = match Sources::read(&root) {
        Ok(s) => s,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            return 1;
        }
    };
    let _ = writeln!(out, "flag\trun_by\tsteps\timages\tsources\tinstalled");
    for c in components::TABLE {
        let sources: Vec<&str> = c.sources(&s).iter().map(|r| r.id.as_str()).collect();
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            c.flag,
            if c.ported() { "og" } else { "bash" },
            c.steps,
            c.images.join(","),
            sources.join(","),
            c.installed.join(" ")
        );
    }
    0
}
