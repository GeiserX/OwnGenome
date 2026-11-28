//! `og setup check [SAMPLE]`: scripts/validate-setup.sh, run by og.
//!
//! It checks what a run needs before it starts, in validate-setup.sh's five
//! sections and with its words: the system (bash, the container engine and
//! its memory, wget or curl, CPUs, RAM, free disk, the architecture), the
//! environment (GENOME_DIR), the reference data, the images of versions.env
//! with their pinned tags, og's own binaries, and with SAMPLE the sample's
//! reads, BAM or VCF. Every line keeps validate-setup.sh's text and its
//! `[OK]`/`[WARN]`/`[FAIL]`/`[INFO]` form, and the exit status is its:
//! 0 with no `[FAIL]`, 1 with one or more, 2 for a sample name or a
//! reference the scripts refuse. run-all.sh, CI and the tests read them.
//!
//! The paths an installer writes come from the component table
//! ([`crate::components`]) that `og setup install` installs from: each
//! check of an installed file asks for it by its row and template
//! ([`Check::table`]), so a row that moves a file moves its check. A path of
//! the table no check reads is listed in [`NOT_CHECKED`] with the reason,
//! and a test fails on a path that is in neither.
//!
//! The sample checks run samtools and bcftools in their images, as
//! validate-setup.sh does (scripts/lib/common.sh `run_in`). Every engine
//! call is the shell's, with the same arguments and in the same order, so a
//! fake engine logs the same lines for both. One difference on purpose: the
//! engine's memory is read from CONTAINER_ENGINE (validate-setup.sh asked a
//! literal `docker`), and from podman's `{{.Host.MemTotal}}` when
//! `{{.MemTotal}}` gives none.

use crate::components;
use crate::engine::{Place, RunIn};
use crate::fetch::{self, Ctx};
use crate::images::{self, Engine};
use crate::install::{Vars, command_v, is_executable, nonempty, sha256};
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::io::{IsTerminal, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

/// The paths of the component table no check reads, each with the reason:
/// (component key, installed template, why). Every other installed path of
/// [`components::TABLE`] is read by a check through [`Check::table`]; a test
/// fails on a path that is in neither, so a new row needs a check or a line
/// here.
pub const NOT_CHECKED: &[(&str, &str, &str)] = &[
    (
        "prs-academic",
        "prs_scores/PGS000014.txt.gz",
        "checked from the academic-research rows of assets/pgs_scores.tsv, the list step 25 reads (a test ties the two)",
    ),
    (
        "prs-academic",
        "prs_scores/PGS000016.txt.gz",
        "checked from the academic-research rows of assets/pgs_scores.tsv (as PGS000014)",
    ),
    (
        "prs-academic",
        "prs_scores/PGS000017.txt.gz",
        "checked from the academic-research rows of assets/pgs_scores.tsv (as PGS000014)",
    ),
    (
        "prs-academic",
        "prs_scores/PGS000027.txt.gz",
        "checked from the academic-research rows of assets/pgs_scores.tsv (as PGS000014)",
    ),
    (
        "parascopy-data",
        "reference/parascopy-${PARASCOPY_DATA_VERSION}/models_GRCh38_1KGP/EUR/SMN1.gz",
        "validate-setup.sh checks the homology table only; the models come with it",
    ),
    (
        "ancestry-fine-panel",
        "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel/DONE",
        "model/DONE, which is checked, is written after it: a finished model means a finished panel",
    ),
    (
        "yleaf-data",
        "reference/yleaf-${YLEAF_DATA_VERSION}/data/hg_prediction_tables/tree.json",
        "validate-setup.sh reads Yleaf's data only through new_positions.txt, to check its image",
    ),
    (
        "ancient-panel",
        "reference/ancient/aadr_${AADR_VERSION}/DONE",
        "opt-in data validate-setup.sh has no line for; main.nf checks this DONE before step 44",
    ),
    (
        "viral-data",
        "reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "masked-reference",
        "{ref_stem}_maskedGRC.fasta",
        "a run reads it only as REF_FASTA, which the reference checks read by that name",
    ),
    (
        "masked-reference",
        "{ref_stem}_maskedGRC.fasta.fai",
        "read as REF_FASTA's .fai when REF_FASTA names the masked reference",
    ),
    (
        "masked-reference",
        "{ref_stem}_maskedGRC.dict",
        "read as REF_FASTA's .dict when REF_FASTA names the masked reference",
    ),
    (
        "exomiser-data",
        "reference/exomiser-${EXOMISER_DATA_VERSION}/${EXOMISER_DATA_VERSION}_hg38",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "exomiser-data",
        "reference/exomiser-${EXOMISER_DATA_VERSION}/${EXOMISER_DATA_VERSION}_phenotype",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "chip-data",
        "reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "strchive-data",
        "reference/strchive-${STRCHIVE_VERSION}/catalog.json",
        "validate-setup.sh has no line for the repeat catalog; a full setup builds it and notes it when missing",
    ),
    (
        "strchive-data",
        "reference/strchive-${STRCHIVE_VERSION}/str_ranges.tsv",
        "as catalog.json",
    ),
    (
        "strchive-data",
        "reference/strchive-${STRCHIVE_VERSION}/catalog.stamp",
        "as catalog.json",
    ),
    (
        "megane-kmers",
        "reference/megane/{ref_name}.mk",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "megane-kmers",
        "reference/megane/{ref_name}.mi",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "mocha-data",
        "reference/mocha-${MOCHA_PANEL_VERSION}/DONE",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "locityper-data",
        "reference/locityper-${LOCITYPER_PANEL_VERSION}/DONE",
        "opt-in data validate-setup.sh has no line for; main.nf checks this DONE before step 53",
    ),
    (
        "locityper-data",
        "reference/locityper-${LOCITYPER_PANEL_VERSION}/counts.jf",
        "as DONE",
    ),
    (
        "locityper-data",
        "reference/locityper-${LOCITYPER_PANEL_VERSION}/loci.bed",
        "as DONE",
    ),
    (
        "locityper-data",
        "reference/locityper-${LOCITYPER_PANEL_VERSION}/db",
        "as DONE",
    ),
    (
        "annotation-scores",
        "annotation_scores",
        "the row names annotation_scores, but setup.sh's install_annotation_scores writes annotations/, which steps 30 and 31 read and the annotation check reads; the annotation port corrects the row",
    ),
    (
        "somatic-data",
        "somatic/af-only-gnomad.hg38.vcf.gz",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "somatic-data",
        "somatic/af-only-gnomad.hg38.vcf.gz.tbi",
        "as af-only-gnomad.hg38.vcf.gz",
    ),
    (
        "somatic-data",
        "somatic/1000g_pon.hg38.vcf.gz",
        "as af-only-gnomad.hg38.vcf.gz",
    ),
    (
        "somatic-data",
        "somatic/1000g_pon.hg38.vcf.gz.tbi",
        "as af-only-gnomad.hg38.vcf.gz",
    ),
    (
        "somatic-data",
        "somatic/small_exac_common_3.hg38.vcf.gz",
        "as af-only-gnomad.hg38.vcf.gz",
    ),
    (
        "somatic-data",
        "somatic/small_exac_common_3.hg38.vcf.gz.tbi",
        "as af-only-gnomad.hg38.vcf.gz",
    ),
    (
        "bwa-index",
        "{ref}.amb",
        "opt-in index validate-setup.sh has no line for",
    ),
    ("bwa-index", "{ref}.ann", "as {ref}.amb"),
    ("bwa-index", "{ref}.bwt", "as {ref}.amb"),
    ("bwa-index", "{ref}.pac", "as {ref}.amb"),
    ("bwa-index", "{ref}.sa", "as {ref}.amb"),
    (
        "cnvpytor-data",
        "reference/cnvpytor/gc_hg38.pytor",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "cnvpytor-data",
        "reference/cnvpytor/mask_hg38.pytor",
        "as gc_hg38.pytor",
    ),
    (
        "cnvpytor-data",
        "reference/cnvpytor/gc_hg19.pytor",
        "as gc_hg38.pytor",
    ),
    (
        "cnvpytor-data",
        "reference/cnvpytor/mask_hg19.pytor",
        "as gc_hg38.pytor",
    ),
    (
        "cnvpytor-data",
        "reference/cnvpytor/gc_chm13v2.0.pytor",
        "as gc_hg38.pytor",
    ),
    (
        "cnvpytor-data",
        "reference/cnvpytor/gc_chm13v1.1.pytor",
        "as gc_hg38.pytor",
    ),
    (
        "cnvpytor-data",
        "reference/cnvpytor/gc_kn99.pytor",
        "as gc_hg38.pytor",
    ),
    (
        "pypgx-data",
        "reference/pypgx-bundle",
        "opt-in data validate-setup.sh has no line for",
    ),
    (
        "archaic",
        "reference/archaic/${ARCHAIC_DATA_VERSION}/DONE",
        "opt-in data validate-setup.sh has no line for; main.nf checks this DONE before step 43",
    ),
];

/// common.sh's DATA_FILES, in its order (a test ties the two).
pub const DATA_FILES: [&str; 5] = [
    "delly_exclude",
    "cytoband",
    "hla_dat",
    "gencode_genes",
    "gencode_exons",
];

/// The colour codes of validate-setup.sh: tput's, only when stdout is a
/// terminal with 8 colours or more.
#[derive(Debug, Clone, Default)]
pub struct Colors {
    pub green: String,
    pub yellow: String,
    pub red: String,
    pub bold: String,
    pub reset: String,
}

impl Colors {
    /// `[ -t 1 ] && command -v tput && [ "$(tput colors)" -ge 8 ]`.
    pub fn for_stdout() -> Colors {
        if !std::io::stdout().is_terminal() || !command_v("tput") {
            return Colors::default();
        }
        let tput = |a: &[&str]| -> (bool, String) { capture("tput", a) };
        let (ok, n) = tput(&["colors"]);
        let n = if ok { n } else { or_zero(n) };
        if test_int(&n).is_none_or(|n| n < 8) {
            return Colors::default();
        }
        Colors {
            green: tput(&["setaf", "2"]).1,
            yellow: tput(&["setaf", "3"]).1,
            red: tput(&["setaf", "1"]).1,
            bold: tput(&["bold"]).1,
            reset: tput(&["sgr0"]).1,
        }
    }
}

/// The report: validate-setup.sh's pass, warn, fail, info and header lines,
/// and its two counters.
pub struct Report<'a> {
    out: &'a mut dyn Write,
    pub c: Colors,
    pub failures: u32,
    pub warnings: u32,
}

impl<'a> Report<'a> {
    pub fn new(out: &'a mut dyn Write, c: Colors) -> Report<'a> {
        Report {
            out,
            c,
            failures: 0,
            warnings: 0,
        }
    }

    /// `echo "LINE"`.
    pub fn say(&mut self, line: &str) {
        let _ = writeln!(self.out, "{line}");
    }

    pub fn pass(&mut self, m: &str) {
        let l = format!("  {}[OK]{}    {m}", self.c.green, self.c.reset);
        self.say(&l);
    }

    pub fn warn(&mut self, m: &str) {
        let l = format!("  {}[WARN]{}  {m}", self.c.yellow, self.c.reset);
        self.say(&l);
        self.warnings += 1;
    }

    pub fn fail(&mut self, m: &str) {
        let l = format!("  {}[FAIL]{}  {m}", self.c.red, self.c.reset);
        self.say(&l);
        self.failures += 1;
    }

    pub fn info(&mut self, m: &str) {
        let l = format!("  {}[INFO]{}  {m}", self.c.bold, self.c.reset);
        self.say(&l);
    }

    pub fn header(&mut self, t: &str) {
        self.say("");
        let l = format!("{}=== {t} ==={}", self.c.bold, self.c.reset);
        self.say(&l);
    }

    pub fn flush(&mut self) {
        let _ = self.out.flush();
    }
}

/// What check_bam_reference found of the BAM's @SQ M5s (BAM_M5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum M5 {
    /// Every @SQ line has an M5 and each is REF_DICT's.
    Match,
    /// One differs.
    Differ,
    /// Anything else.
    Unknown,
}

/// Everything one check reads: the checkout, the environment as common.sh
/// leaves it, and versions.env.
pub struct Check {
    pub root: PathBuf,
    pub vars: Vars,
    pub ctx: Ctx,
    pub engine: Engine,
    /// GENOME_DIR without its trailing slash; "" when unset.
    pub genome: String,
    pub sample: String,
    pub ref_fasta: String,
    pub ref_dict: String,
    pub ref_fasta_c: String,
    /// The (component, template) pairs a check asked the table for.
    pub seen: RefCell<BTreeSet<(String, String)>>,
}

/// `[ -f P ]`.
fn is_file(p: &str) -> bool {
    std::fs::metadata(p).is_ok_and(|m| m.is_file())
}

/// `[ -d P ]`.
fn is_dir(p: &str) -> bool {
    std::fs::metadata(p).is_ok_and(|m| m.is_dir())
}

/// `[ -e P ]`.
fn exists(p: &str) -> bool {
    std::fs::metadata(p).is_ok()
}

/// `[ -s P ]`.
fn s(p: &str) -> bool {
    nonempty(Path::new(p))
}

fn mtime(p: &str) -> Option<SystemTime> {
    std::fs::metadata(p).ok().and_then(|m| m.modified().ok())
}

/// `[ A -ot B ]`: A is older than B, or B exists and A does not.
pub fn older(a: &str, b: &str) -> bool {
    match (mtime(a), mtime(b)) {
        (Some(x), Some(y)) => x < y,
        (None, Some(_)) => true,
        _ => false,
    }
}

/// `[ A -nt B ]`: A is newer than B, or A exists and B does not.
pub fn newer(a: &str, b: &str) -> bool {
    match (mtime(a), mtime(b)) {
        (Some(x), Some(y)) => x > y,
        (Some(_), None) => true,
        _ => false,
    }
}

/// The lines of TEXT as awk and grep read them: a last line without its
/// newline counts, the empty string after a final newline does not.
pub fn lines(text: &str) -> Vec<&str> {
    if text.is_empty() {
        return vec![];
    }
    let t = text.strip_suffix('\n').unwrap_or(text);
    t.split('\n').collect()
}

/// `grep -c .`: the lines with at least one character.
pub fn count_nonempty(text: &str) -> usize {
    lines(text).iter().filter(|l| !l.is_empty()).count()
}

/// `$(cat FILE 2>/dev/null)`: the content without its trailing newlines,
/// "" when it cannot be read.
fn cat(p: &str) -> String {
    std::fs::read(p)
        .map(|b| {
            String::from_utf8_lossy(&b)
                .trim_end_matches('\n')
                .to_owned()
        })
        .unwrap_or_default()
}

fn read_text(p: &str) -> String {
    std::fs::read(p)
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .unwrap_or_default()
}

/// `[ "$X" -gt 0 ]` and friends read X as bash's test does: blanks around
/// an optional sign and decimal digits.
pub fn test_int(x: &str) -> Option<i64> {
    let t = x.trim_matches(|c: char| c == ' ' || c == '\t' || c == '\n');
    let (neg, d) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if d.is_empty() || !d.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let n: i64 = d.parse().ok()?;
    Some(if neg { -n } else { n })
}

/// `[[ X =~ ^[0-9]+$ ]]`.
fn digits(x: &str) -> bool {
    !x.is_empty() && x.bytes().all(|b| b.is_ascii_digit())
}

/// What `$(CMD || echo 0)` gives when CMD failed: its output, then 0.
fn or_zero(out: String) -> String {
    if out.is_empty() {
        "0".to_owned()
    } else {
        format!("{out}\n0")
    }
}

/// `$(NAME ARGS 2>/dev/null)`: stdout without its trailing newlines, and
/// whether the command ran and succeeded.
pub fn capture(name: &str, args: &[&str]) -> (bool, String) {
    match Command::new(name)
        .args(args)
        .stdin(Stdio::null())
        .stderr(Stdio::null())
        .output()
    {
        Ok(o) => (
            o.status.success(),
            String::from_utf8_lossy(&o.stdout)
                .trim_end_matches('\n')
                .to_owned(),
        ),
        Err(_) => (false, String::new()),
    }
}

/// awk's `printf "%.Nf"`.
pub fn fixed(x: f64, n: usize) -> String {
    format!("{x:.n$}")
}

/// `dirname P`.
pub fn dirname(p: &str) -> String {
    let t = p.trim_end_matches('/');
    if t.is_empty() {
        return if p.starts_with('/') { "/" } else { "." }.to_owned();
    }
    match t.rfind('/') {
        None => ".".to_owned(),
        Some(i) => {
            let d = t[..i].trim_end_matches('/');
            if d.is_empty() { "/" } else { d }.to_owned()
        }
    }
}

/// `basename P`.
fn basename(p: &str) -> &str {
    let t = p.trim_end_matches('/');
    t.rsplit('/').next().unwrap_or(t)
}

/// `cd -- DIR && pwd -P`: the physical path, "" when DIR cannot be entered.
fn physical(d: &str) -> String {
    if d.is_empty() {
        return String::new();
    }
    std::fs::canonicalize(d)
        .ok()
        .filter(|p| p.is_dir())
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// fastp_total_reads JSON of common.sh: summary.after_filtering.total_reads,
/// the first "after_filtering" object of the report with every blank
/// removed; `None` when the report has none.
pub fn fastp_total_reads(json: &str) -> Option<String> {
    let b = std::fs::read(json).ok()?;
    let t: String = String::from_utf8_lossy(&b)
        .chars()
        .filter(|c| !matches!(c, ' ' | '\t' | '\r' | '\n'))
        .collect();
    let k = "\"after_filtering\":{\"total_reads\":";
    let i = t.find(k)?;
    let v: String = t[i + k.len()..]
        .chars()
        .take_while(char::is_ascii_digit)
        .collect();
    (!v.is_empty()).then_some(v)
}

/// read_count_check PRIMARY SOURCE "" COUNT of common.sh, its standard
/// output and error in one text without the trailing newline, and its exit
/// status: 0 equal, 1 different, 3 different under READ_COUNT_CHECK=warn, 2
/// for a count that is not a whole number or a bad READ_COUNT_CHECK.
pub fn read_count_check(
    primary: &str,
    src: &str,
    counts: &[&str],
    read_count_check: Option<&str>,
) -> (u8, String) {
    let mut msg: Vec<String> = Vec::new();
    let mode = match read_count_check.filter(|v| !v.is_empty()).unwrap_or("fail") {
        m @ ("fail" | "warn") => m.to_owned(),
        _ => {
            msg.push(format!(
                "ERROR: READ_COUNT_CHECK must be 'fail' or 'warn', got '{}'",
                read_count_check.unwrap_or("")
            ));
            return (2, msg.join("\n"));
        }
    };
    let whole =
        |c: &str| c == "0" || (c.starts_with(|x: char| ('1'..='9').contains(&x)) && digits(c));
    if !std::iter::once(primary)
        .chain(counts.iter().copied())
        .all(whole)
    {
        msg.push(format!(
            "ERROR: could not count the reads: '{primary}' primary records in the BAM, '{}' reads in ({src})",
            counts.join(" ")
        ));
        return (2, msg.join("\n"));
    }
    if counts.is_empty() {
        return (2, "ERROR: read_count_check got no input count".to_owned());
    }
    let wrap = |s: &str| s.parse::<i64>().unwrap_or(0);
    let total: i64 = counts.iter().map(|c| wrap(c)).fold(0i64, i64::wrapping_add);
    let p = wrap(primary);
    let missing = total.wrapping_sub(p);
    if missing == 0 {
        return (
            0,
            format!(
                "Read count: the BAM holds {p} primary records, the {total} reads that went in ({src})."
            ),
        );
    }
    let level = if mode == "fail" { "ERROR" } else { "WARNING" };
    if missing > 0 {
        let pct = fixed(100.0 * missing as f64 / total as f64, 1);
        msg.push(format!(
            "{level}: the BAM holds {p} primary records, but {total} reads went in ({src}): {missing} reads ({pct}%) are missing."
        ));
    } else {
        msg.push(format!(
            "{level}: the BAM holds {p} primary records, more than the {total} reads that went in ({src}): it was not made from these reads."
        ));
    }
    msg.push("  The aligner writes every read it is given, mapped or not, once, so the two counts must be equal.".to_owned());
    msg.push("  Depth, copy number (CYP2D6), CNV, telomere content and every call that depends on coverage would be off.".to_owned());
    msg.push("  An R1 and R2 with different numbers of records is one cause: minimap2 and bwa-mem2 stop at the end".to_owned());
    msg.push("  of the shorter file with only a warning. See docs/troubleshooting.md, \"The BAM lost reads\".".to_owned());
    if mode == "warn" {
        msg.push("  READ_COUNT_CHECK=warn: going on.".to_owned());
        return (3, msg.join("\n"));
    }
    msg.push("  For an input you know is odd, READ_COUNT_CHECK=warn goes on anyway.".to_owned());
    (1, msg.join("\n"))
}

/// vcf_header_samples of common.sh: the columns after FORMAT of every
/// #CHROM line, one per line, as `$(...)` keeps them (trailing newlines
/// dropped).
pub fn vcf_header_samples(header: &str) -> String {
    let mut v: Vec<&str> = Vec::new();
    for l in lines(header) {
        if l.starts_with("#CHROM") {
            v.extend(l.split('\t').skip(9));
        }
    }
    v.join("\n").trim_end_matches('\n').to_owned()
}

/// validate_sample of common.sh: the message, `Err` when SAMPLE is refused.
pub fn validate_sample(s: &str) -> Result<(), String> {
    let plain = !s.is_empty()
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-'));
    if !plain || s == "." || s == ".." {
        return Err(format!(
            "ERROR: invalid sample name '{s}'. Use only letters, digits, '.', '_' and '-' (and not '.' or '..')."
        ));
    }
    if s.starts_with('-') {
        return Err(format!(
            "ERROR: invalid sample name '{s}': it starts with '-', which tools read as an option. Start it with a letter or digit."
        ));
    }
    Ok(())
}

/// The fields of a SAM header line after its record type: TAG -> value of
/// the first field with that tag (awk's loop that keeps the last one is
/// [`last_tag`]).
fn first_tag<'a>(fields: &[&'a str], t: &str) -> &'a str {
    let k = format!("{t}:");
    fields
        .iter()
        .skip(1)
        .find_map(|f| f.strip_prefix(k.as_str()))
        .unwrap_or("")
}

/// awk's `for (i = 2; i <= NF; i++) if ($i ~ /^T:/) v = substr($i, 4)`: the
/// last field with the two-letter tag T.
fn last_tag<'a>(fields: &[&'a str], t: &str) -> &'a str {
    let k = format!("{t}:");
    fields
        .iter()
        .skip(1)
        .filter_map(|f| f.strip_prefix(k.as_str()))
        .next_back()
        .unwrap_or("")
}

impl Check {
    /// What common.sh sets up once sourced, from VARS (the process
    /// environment): GENOME_DIR without its trailing slash, the downloads
    /// folder, REF_FASTA (absolute, inside GENOME_DIR), REF_DICT and
    /// REF_FASTA_C, CONTAINER_ENGINE and THREADS; SAMPLE is the argument.
    /// `Err` carries the exit code and the message of a reference outside
    /// GENOME_DIR.
    pub fn new(root: &Path, mut vars: Vars, sample: &str) -> Result<Check, (u8, String)> {
        let genome_dir = vars.get("GENOME_DIR").unwrap_or("").to_owned();
        let mut msg = Vec::new();
        let run = crate::install::Run::new(root, vars.clone(), "", &genome_dir, true, &mut msg)
            .map_err(|c| (c, String::from_utf8_lossy(&msg).into_owned()))?;
        vars = run.vars;
        vars.set("SAMPLE", sample);
        let genome = run.genome;
        let (ref_fasta, ref_dict, ref_fasta_c) = if genome.is_empty() {
            (String::new(), String::new(), String::new())
        } else {
            let r = vars.get("REF_FASTA").unwrap_or("").to_owned();
            let d = r.strip_suffix(".gz").unwrap_or(&r);
            let d = match d.rfind('.') {
                Some(i) => format!("{}.dict", &d[..i]),
                None => format!("{d}.dict"),
            };
            let c = crate::engine::cpath(&genome, vars.nonempty("OG_DOWNLOADS"), &r)
                .map_err(|m| (2u8, format!("{m}\n")))?;
            (r, d, c)
        };
        let ctx = Ctx::new(root, vars.boxed()).map_err(|m| (1u8, format!("{m}\n")))?;
        let engine = Engine::from_vars(&vars);
        Ok(Check {
            root: root.to_path_buf(),
            vars,
            ctx,
            engine,
            genome,
            sample: sample.to_owned(),
            ref_fasta,
            ref_dict,
            ref_fasta_c,
            seen: RefCell::new(BTreeSet::new()),
        })
    }

    /// A versions.env value, or a setting of common.sh (the environment's
    /// when common.sh takes it from there); "" when neither has it.
    pub fn v(&self, name: &str) -> String {
        self.ctx
            .expand(&format!("${{{name}}}"), None)
            .unwrap_or_default()
    }

    /// An environment variable, "" when unset.
    fn env(&self, name: &str) -> &str {
        self.vars.get(name).unwrap_or("")
    }

    /// The installed path TEMPLATE of the component KEY, expanded: the one
    /// way a check names a file an installer writes. Panics when the row
    /// has no such template (a check that drifted from the table, caught by
    /// the tests).
    pub fn table(&self, key: &str, template: &str) -> String {
        let row = components::by_key(key)
            .unwrap_or_else(|| panic!("og setup check: no component {key} in the table"));
        let i = row
            .installed
            .iter()
            .position(|t| *t == template)
            .unwrap_or_else(|| {
                panic!("og setup check: the {key} row has no installed path {template}")
            });
        self.seen
            .borrow_mut()
            .insert((key.to_owned(), template.to_owned()));
        let paths = row
            .installed_paths(
                &self.ctx,
                Path::new(&self.genome),
                Path::new(&self.ref_fasta),
            )
            .unwrap_or_default();
        paths
            .get(i)
            .map(|p| p.to_string_lossy().into_owned())
            .unwrap_or_else(|| format!("{}/{template}", self.genome))
    }

    /// The steps of KEY's row.
    fn steps(key: &str) -> &'static str {
        components::by_key(key).map_or("", |c| c.steps)
    }

    /// The flag of KEY's row.
    fn flag(key: &str) -> &'static str {
        components::by_key(key).map_or("", |c| c.flag)
    }

    /// The image variable VAR of KEY's row (a check that drifted from the
    /// table panics, as in [`Check::table`]).
    fn image_var(key: &str, var: &'static str) -> &'static str {
        let row = components::by_key(key)
            .unwrap_or_else(|| panic!("og setup check: no component {key} in the table"));
        assert!(
            row.images.contains(&var),
            "og setup check: the {key} row pulls no {var}"
        );
        var
    }

    /// run_in IMAGE COMMAND... with its output kept and its errors dropped:
    /// the exit status (1 when it could not start) and stdout without its
    /// trailing newlines.
    fn run_in(&self, rest: &[&str]) -> (i32, String) {
        let place = Place::from_vars(&self.vars);
        let rest: Vec<String> = rest.iter().map(|s| (*s).to_owned()).collect();
        let Ok(a) = RunIn::default().argv(&place, &rest) else {
            return (1, String::new());
        };
        match Command::new(&self.engine.cmd)
            .args(&a)
            .stdin(Stdio::null())
            .stderr(Stdio::null())
            .output()
        {
            Ok(o) => (
                o.status.code().unwrap_or(1),
                String::from_utf8_lossy(&o.stdout)
                    .trim_end_matches('\n')
                    .to_owned(),
            ),
            Err(_) => (1, String::new()),
        }
    }

    fn bam_c(&self) -> String {
        format!("/genome/{0}/aligned/{0}_sorted.bam", self.sample)
    }

    /// check_bam_reference HEADER: the BAM's @SQ names and lengths are the
    /// .fai's in the same order, its M5s REF_DICT's, and it is
    /// coordinate-sorted.
    pub fn check_bam_reference(&self, r: &mut Report, header: &str) -> M5 {
        let fai = format!("{}.fai", self.ref_fasta);
        let mut bam_m5 = M5::Unknown;
        let hl = lines(header);
        if is_file(&fai) && hl.iter().any(|l| l.starts_with("@SQ")) {
            let fai_text = read_text(&fai);
            let bam: Vec<String> = hl
                .iter()
                .filter(|l| l.starts_with("@SQ"))
                .map(|l| {
                    let f: Vec<&str> = l.split('\t').collect();
                    format!("{} ({} bp)", last_tag(&f, "SN"), last_tag(&f, "LN"))
                })
                .collect();
            let reference: Vec<String> = lines(&fai_text)
                .iter()
                .filter_map(|l| {
                    let f: Vec<&str> = l.split('\t').collect();
                    (f.len() >= 2).then(|| format!("{} ({} bp)", f[0], f[1]))
                })
                .collect();
            let n = bam.len().max(reference.len());
            let mut diff = String::new();
            for i in 0..n {
                if bam.get(i) != reference.get(i) {
                    diff = format!(
                        "sequence {} is {} in the BAM and {} in the reference",
                        i + 1,
                        bam.get(i).map_or("missing", String::as_str),
                        reference.get(i).map_or("missing", String::as_str)
                    );
                    break;
                }
            }
            if diff.is_empty() {
                r.pass(&format!(
                    "BAM header matches the reference: the same {} sequences in the same order",
                    count_nonempty(&fai_text)
                ));
                if s(&self.ref_dict) {
                    match self.m5_compare(&hl) {
                        Some(Ok(())) => {
                            bam_m5 = M5::Match;
                            r.pass("BAM's @SQ MD5s (M5) match the reference dictionary: it was aligned to these bases");
                        }
                        Some(Err(d)) => {
                            bam_m5 = M5::Differ;
                            r.fail("this BAM was aligned to a FASTA with other bases than REF_FASTA, though with the same names and lengths (the masked reference, or another copy of GRCh38): realign (docs/realignment.md)");
                            r.say(&format!("       First difference: {d}."));
                            r.say(&format!("       Reference: {}", self.ref_fasta));
                        }
                        None => {}
                    }
                }
            } else {
                r.fail(
                    "this BAM was aligned to a different reference: realign (docs/realignment.md)",
                );
                r.say(&format!("       First difference: {diff}."));
                r.say(&format!("       Reference: {}", self.ref_fasta));
            }
        }
        let so: Vec<&str> = hl
            .iter()
            .filter(|l| l.starts_with("@HD"))
            .flat_map(|l| l.split('\t').skip(1).filter_map(|f| f.strip_prefix("SO:")))
            .collect();
        let so = so.join("\n");
        match so.as_str() {
            "coordinate" => r.pass("BAM is coordinate-sorted (@HD SO:coordinate)"),
            "" => {
                let bam = format!(
                    "{}/{}/aligned/{}_sorted.bam",
                    self.genome, self.sample, self.sample
                );
                let bai = format!("{bam}.bai");
                if is_file(&bai) && !older(&bai, &bam) {
                    r.pass("BAM is coordinate-sorted (no @HD SO: tag, but samtools only indexes a sorted BAM and its .bai is not older than it)");
                } else if is_file(&bai) {
                    r.fail("BAM header does not say how it is sorted (no @HD SO: tag) and its .bai is older than the BAM: index it again (samtools index fails on an unsorted BAM) or realign with step 02");
                } else {
                    r.fail("BAM header does not say how it is sorted (no @HD SO: tag) and it has no .bai: sort and index it (samtools sort, samtools index) or realign with step 02");
                }
            }
            other => r.fail(&format!(
                "BAM is sorted by {other}, not by coordinate: sort it (samtools sort) or realign with step 02"
            )),
        }
        bam_m5
    }

    /// The M5 half of check_bam_reference: `Some(Ok)` when every @SQ line
    /// of the BAM has an M5 and each is REF_DICT's, `Some(Err(first
    /// difference))` on a difference, `None` otherwise.
    fn m5_compare(&self, header: &[&str]) -> Option<Result<(), String>> {
        let dict = read_text(&self.ref_dict);
        let mut reference: std::collections::HashMap<String, String> =
            std::collections::HashMap::new();
        for l in lines(&dict) {
            if l.starts_with("@SQ") {
                let f: Vec<&str> = l.split('\t').collect();
                let (sn, m5) = (first_tag(&f, "SN"), first_tag(&f, "M5").to_lowercase());
                if !sn.is_empty() && !m5.is_empty() {
                    reference.insert(sn.to_owned(), m5);
                }
            }
        }
        let (mut n, mut hit, mut partial) = (0usize, 0usize, false);
        for l in header {
            if !l.starts_with("@SQ") {
                continue;
            }
            n += 1;
            let f: Vec<&str> = l.split('\t').collect();
            let (sn, b) = (first_tag(&f, "SN"), first_tag(&f, "M5").to_lowercase());
            match reference.get(sn) {
                Some(want) if !b.is_empty() => {
                    if &b != want {
                        return Some(Err(format!(
                            "sequence {n} ({sn}) has MD5 {b} in the BAM and {want} in the reference dictionary"
                        )));
                    }
                    hit += 1;
                }
                _ => partial = true,
            }
        }
        (hit > 0 && !partial).then_some(Ok(()))
    }

    /// check_bam_masked_reference HEADER: the BAM was aligned to the masked
    /// reference exactly when REF_FASTA is that reference, as its aligner's
    /// @PG lines name it.
    pub fn check_bam_masked_reference(&self, r: &mut Report, header: &str, bam_m5: M5) {
        let want = basename(&self.ref_fasta).contains("_maskedGRC.");
        let (mut default, mut masked) = (false, false);
        for l in lines(header) {
            if !l.starts_with("@PG") {
                continue;
            }
            let f: Vec<&str> = l.split('\t').collect();
            let mut id = last_tag(&f, "ID").to_owned();
            let pn = last_tag(&f, "PN");
            let cl = last_tag(&f, "CL");
            // samtools merge writes ID:minimap2.1, a second merge minimap2.1.1.
            while let Some(i) = id.rfind('.')
                && digits(&id[i + 1..])
            {
                id.truncate(i);
            }
            let aligner = |x: &str| matches!(x, "minimap2" | "bwa-mem2" | "bwa");
            if !aligner(&id) && !aligner(pn) {
                continue;
            }
            if cl.is_empty() {
                continue;
            }
            let mut kind_masked = false;
            for w in cl.split(' ') {
                let x = w.trim_matches(|c| c == '"' || c == '\'');
                if x.starts_with(['-', '@']) || read_file_word(x) {
                    continue;
                }
                if basename_any(x).contains("_maskedGRC") {
                    kind_masked = true;
                }
            }
            if kind_masked {
                masked = true;
            } else {
                default = true;
            }
        }
        let both = default && masked;
        if !both && bam_m5 == M5::Match && !older(&self.ref_dict, &self.ref_fasta) {
            return;
        }
        match (want, default, masked) {
            (_, true, true) => r.fail("this BAM's aligner @PG lines name both the masked and the default reference (reads merged from the two): realign it to one (docs/realignment.md#the-masked-reference)"),
            (true, false, false) => {
                r.warn("BAM header has no minimap2 or BWA @PG line with a command line, so it cannot show it was aligned to the masked reference REF_FASTA names");
                r.say("       If it was aligned to the default reference, realign it (docs/realignment.md#the-masked-reference).");
            }
            (true, true, false) => r.fail("REF_FASTA is the masked reference, but this BAM was aligned to the default one (its aligner @PG line names no _maskedGRC file): realign it (docs/realignment.md#the-masked-reference)"),
            (false, false, true) => r.fail(&format!(
                "this BAM was aligned to the masked reference (its aligner @PG line names a _maskedGRC file), but REF_FASTA is {}: set REF_FASTA to the masked reference (docs/realignment.md#the-masked-reference)",
                self.ref_fasta
            )),
            (true, false, true) => r.pass("BAM was aligned to the masked reference (its aligner @PG line names a _maskedGRC file)"),
            _ => {}
        }
    }

    /// check_reference_content: REF_FASTA holds the bases of NCBI's GRCh38
    /// no-ALT analysis set, or of its masked copy when its name says so, by
    /// the M5s REF_DICT holds.
    pub fn check_reference_content(&self, r: &mut Report) {
        const SHAPE: &str = "8502b755744bd2d708d866bf3e68a58c7bd9369c737e922ba32d7fea4aedaf51";
        const STOCK: &str = "f6d0a052602c69fe61c3efe9aab12112f1cafc2b4a5140d4fd240d27e6ba8743";
        const STOCK_MASKED: &str =
            "73724d8ce85ff35460f369c9d0cd899a285258e1e66dcbaeedf22ab90518712c";
        let dict = &self.ref_dict;
        let fai = format!("{}.fai", self.ref_fasta);
        if !(s(dict) && is_file(&fai)) {
            return;
        }
        let g = &self.genome;
        let masked = basename(&self.ref_fasta).contains("_maskedGRC.");
        let rows: Vec<(String, String, String)> = lines(&read_text(dict))
            .iter()
            .filter(|l| l.starts_with("@SQ"))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (
                    last_tag(&f, "SN").to_owned(),
                    last_tag(&f, "LN").to_owned(),
                    last_tag(&f, "M5").to_lowercase(),
                )
            })
            .collect();
        // $(cut -f1,2 <<< "$lines"): a here-string of no rows is one empty line.
        let shape_lines: Vec<String> = if rows.is_empty() {
            vec![String::new()]
        } else {
            rows.iter().map(|(n, l, _)| format!("{n}\t{l}")).collect()
        };
        let fai_text = read_text(&fai);
        let fai_shape: Vec<String> = lines(&fai_text)
            .iter()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                if f.len() > 1 {
                    format!("{}\t{}", f[0], f[1])
                } else {
                    (*l).to_owned()
                }
            })
            .collect();
        let trim = |v: &[String]| v.join("\n").trim_end_matches('\n').to_owned();
        if trim(&shape_lines) != trim(&fai_shape) {
            r.fail(&format!(
                "The sequence dictionary {dict} does not list the names and lengths of {fai}: it was made from another FASTA"
            ));
            r.say(&format!(
                "       Delete it and run ./scripts/setup.sh {g}, which makes it again from REF_FASTA."
            ));
            return;
        }
        if older(dict, &self.ref_fasta) {
            r.warn(&format!(
                "The sequence dictionary {dict} is older than REF_FASTA, so the MD5s it holds may be of an earlier FASTA"
            ));
            r.say(&format!(
                "       If the FASTA changed, delete the .dict and run ./scripts/setup.sh {g}; if it was only copied, touch the .dict."
            ));
        }
        let digest = |ls: &[String]| {
            let mut h = <sha2::Sha256 as sha2::Digest>::new();
            for l in ls {
                sha2::Digest::update(&mut h, l.as_bytes());
                sha2::Digest::update(&mut h, b"\n");
            }
            hex(&sha2::Digest::finalize(h))
        };
        if digest(&shape_lines) != SHAPE {
            r.info("Reference is not NCBI's GRCh38 no-ALT analysis set (other sequence names or lengths): no pinned MD5s to check its bases against");
            return;
        }
        if rows.iter().any(|(_, _, m)| m.is_empty()) {
            r.warn(&format!(
                "The sequence dictionary {dict} has no MD5 (M5) for some sequences, so the reference's bases are not checked"
            ));
            r.say(&format!(
                "       Delete it and run ./scripts/setup.sh {g}, which makes it again with them."
            ));
            return;
        }
        let full: Vec<String> = rows
            .iter()
            .map(|(n, l, m)| format!("{n}\t{l}\t{m}"))
            .collect();
        let bases = digest(&full);
        let nseq = count_nonempty(&fai_text);
        match (bases.as_str(), masked) {
            (STOCK, false) => r.pass(&format!(
                "Reference bases: NCBI's GRCh38 no-ALT analysis set (the MD5 of each of its {nseq} sequences in the .dict is the pinned one)"
            )),
            (STOCK_MASKED, true) => r.pass(&format!(
                "Reference bases: the masked copy of NCBI's GRCh38 no-ALT analysis set (the MD5 of each of its {nseq} sequences in the .dict is the pinned one)"
            )),
            (STOCK, true) => r.fail(&format!(
                "REF_FASTA is named as the masked reference (_maskedGRC) but holds the bases of the default one: build it again (./scripts/setup.sh --masked-reference {g} after deleting it, its .fai and its .dict)"
            )),
            (STOCK_MASKED, false) => r.fail("REF_FASTA holds the bases of the masked reference but its name has no _maskedGRC, which the BAM check goes by: use the file setup.sh --masked-reference writes (docs/realignment.md#the-masked-reference)"),
            _ => {
                r.fail(&format!(
                    "REF_FASTA has the names and lengths of NCBI's GRCh38 no-ALT analysis set but other bases: the MD5s in {dict} match neither that set nor its masked copy"
                ));
                if masked {
                    r.say(&format!(
                        "       Build it again: delete {}, its .fai and its .dict, and run ./scripts/setup.sh --masked-reference {g}.",
                        self.ref_fasta
                    ));
                } else {
                    r.say(&format!(
                        "       Download it again: delete {}, its .fai and its .dict, and run ./scripts/setup.sh {g}.",
                        self.ref_fasta
                    ));
                }
            }
        }
    }

    /// check_stale_indexes DIR: a warning for each index (.bai, .crai, .csi,
    /// .tbi) under DIR, outside DIR/nextflow, older than its data file, in
    /// the order `find` walks the folder.
    pub fn check_stale_indexes(&self, r: &mut Report, dir: &str) {
        let mut found = Vec::new();
        walk(Path::new(dir), &format!("{dir}/nextflow"), &mut found);
        let mut n = 0;
        for idx in found {
            let mut data = match idx.rfind('.') {
                Some(i) if !idx[i..].contains('/') => idx[..i].to_owned(),
                _ => idx.clone(),
            };
            if !exists(&data) {
                if let Some(b) = idx.strip_suffix(".bai") {
                    data = format!("{b}.bam");
                } else if let Some(b) = idx.strip_suffix(".crai") {
                    data = format!("{b}.cram");
                }
            }
            if !(exists(&data) && older(&idx, &data)) {
                continue;
            }
            n += 1;
            let shown = idx
                .strip_prefix(&format!("{}/", self.genome))
                .unwrap_or(&idx);
            r.warn(&format!("Index older than its data file: {shown}"));
        }
        if n > 0 {
            r.say("       If the data file changed after its index was made, index it again (samtools index for a BAM or CRAM, bcftools index for a VCF); if it was only copied, touch the index.");
        }
    }

    /// check_bam_quickcheck.
    fn check_bam_quickcheck(&self, r: &mut Report) {
        let img = self.v("SAMTOOLS_IMAGE");
        let bam = self.bam_c();
        if self.run_in(&[&img, "samtools", "quickcheck", &bam]).0 == 0 {
            r.pass("BAM passes samtools quickcheck");
        } else {
            r.fail("BAM fails samtools quickcheck (truncated, no sequences in its header, or not a BAM): realign with step 02");
        }
    }

    /// check_bam_read_count: the BAM's primary records equal the reads fastp
    /// wrote, when fastp's report is there.
    fn check_bam_read_count(&self, r: &mut Report, sample_dir: &str) {
        let sm = &self.sample;
        let json = format!("{sample_dir}/fastq_trimmed/{sm}_fastp.json");
        let table = format!("{sample_dir}/qc/{sm}_read_count.tsv");
        let bam = format!("{sample_dir}/aligned/{sm}_sorted.bam");
        if !is_file(&json) {
            return;
        }
        if is_file(&table) && !older(&table, &bam) {
            let src = lines(&read_text(&table))
                .iter()
                .find_map(|l| {
                    let mut f = l.split('\t');
                    (f.next() == Some("input_source")).then(|| f.next().unwrap_or("").to_owned())
                })
                .unwrap_or_default();
            if src.starts_with("fastq/") {
                r.info(&format!(
                    "BAM was aligned from the raw reads (qc/{sm}_read_count.tsv: {src}), not from fastp's output: its read count is not compared with fastp's"
                ));
                return;
            }
        }
        let Some(reads) = fastp_total_reads(&json) else {
            r.warn(&format!(
                "fastq_trimmed/{sm}_fastp.json has no summary.after_filtering.total_reads: the BAM's read count is not checked"
            ));
            return;
        };
        r.info(&format!(
            "Counting the BAM's primary records against fastp's {reads} reads (a few minutes on a 30x BAM)..."
        ));
        r.flush();
        let img = self.v("SAMTOOLS_IMAGE");
        let threads = self.vars.nonempty("THREADS").unwrap_or("8").to_owned();
        let bam_c = self.bam_c();
        let (rc, out) = self.run_in(&[
            "--memory", "2g", &img, "samtools", "view", "-c", "-F", "0x900", "-@", &threads, &bam_c,
        ]);
        let primary = if rc == 0 { out } else { String::new() };
        let src = format!("fastq_trimmed/{sm}_fastp.json");
        let (rc, msg) = read_count_check(
            &primary,
            &src,
            &[reads.as_str()],
            self.vars.get("READ_COUNT_CHECK"),
        );
        match rc {
            0 => r.pass(&format!(
                "BAM holds the {reads} reads fastp wrote as primary records"
            )),
            1 => r.fail(&format!(
                "BAM read count: {primary} primary records, but fastp wrote {reads} reads"
            )),
            3 => r.warn(&format!(
                "BAM read count: {primary} primary records, but fastp wrote {reads} reads (READ_COUNT_CHECK=warn)"
            )),
            _ => r.warn("Could not compare the BAM's read count with fastp's report"),
        }
        if rc != 0 {
            for l in msg.split('\n') {
                r.say(&format!("       {}", l.trim_start_matches(' ')));
            }
            if rc == 1 || rc == 3 {
                r.say("       Realign with step 02 (docs/troubleshooting.md, \"The BAM lost reads\"). A BAM aligned from the");
                r.say(&format!(
                    "       untrimmed reads (FASTQ_SUBDIR=fastq) before step 02 wrote qc/{sm}_read_count.tsv is not"
                ));
                r.say("       comparable with fastp's count: realigning records its own count, and this check then skips it.");
            }
        }
    }

    /// check_vcf_ref_bases: the REF bases of the VCF's first 1000 records
    /// against REF_FASTA (bcftools norm -c x).
    fn check_vcf_ref_bases(&self, r: &mut Report) {
        let v = format!("/genome/{0}/vcf/{0}.vcf.gz", self.sample);
        let script = format!(
            "set -e\n    bcftools view -h '{v}' > /tmp/in.vcf\n    bcftools view -H '{v}' | head -n 1000 >> /tmp/in.vcf\n    bcftools norm -c x -f '{}' -o /tmp/kept.vcf /tmp/in.vcf 2>/dev/null\n    echo \"$(grep -vc '^#' /tmp/in.vcf) $(grep -vc '^#' /tmp/kept.vcf)\"",
            self.ref_fasta_c
        );
        let img = self.v("BCFTOOLS_IMAGE");
        let (rc, out) = self.run_in(&[&img, "sh", "-c", &script]);
        let counts = if rc == 0 { out } else { String::new() };
        let (n, kept) = vcf_counts(&counts);
        r.vcf_ref_bases_verdict(n.as_deref(), kept.as_deref());
    }

    /// The whole check; its exit status.
    pub fn run(&self, r: &mut Report) -> u8 {
        let free_gb = self.system(r);
        self.environment(r);
        self.reference(r);
        let missing_images = self.images(r);
        self.og_binary(r);
        if !self.sample.is_empty() {
            self.sample_data(r, free_gb);
        }
        summary(r, missing_images);
        u8::from(r.failures != 0)
    }

    // --- 1. System Requirements -------------------------------------------------------------

    /// The system section; the free GB of GENOME_DIR's disk, which the
    /// sample section reads.
    fn system(&self, r: &mut Report) -> Option<i64> {
        r.header("System Requirements");
        // 4.4: run-all.sh refuses anything older (it expands empty arrays
        // under set -u). The bash on PATH, which runs run-all.sh.
        let (_, v) = capture(
            "bash",
            &[
                "-c",
                "printf '%s %s' \"${BASH_VERSINFO[0]}\" \"${BASH_VERSINFO[1]}\"",
            ],
        );
        let mut it = v.split(' ');
        let maj = it.next().unwrap_or("").to_owned();
        let min = it.next().unwrap_or("").to_owned();
        let (a, b) = (
            maj.parse::<i64>().unwrap_or(0),
            min.parse::<i64>().unwrap_or(0),
        );
        if a * 100 + b >= 404 {
            r.pass(&format!("bash {maj}.{min} (>= 4.4 required)"));
        } else {
            r.fail(&format!("bash {maj}.{min} — version 4.4+ is required. Install a newer bash (on macOS: brew install bash)."));
        }

        let e = &self.engine;
        let have = command_v(&e.cmd);
        if have {
            let (_, ver) = capture(&e.cmd, &["--version"]);
            let first = ver.lines().next().unwrap_or("").to_owned();
            r.pass(&format!("Docker installed: {first}"));
        } else {
            r.fail("Docker is not installed. Install from https://docs.docker.com/get-docker/");
        }

        let mut engine_mem = String::new();
        if have && e.quiet(&["info"]) {
            r.pass("Docker daemon is running");
            engine_mem = self.engine_memory();
            if test_int(&engine_mem).is_some_and(|n| n > 0) {
                let bytes: f64 = engine_mem.trim().parse().unwrap_or(0.0);
                let gb = fixed(bytes / 1_073_741_824.0, 1);
                let g: f64 = gb.parse().unwrap_or(0.0);
                if g >= 16.0 {
                    r.pass(&format!("Docker memory: {gb} GB (>= 16 GB recommended)"));
                } else if g >= 8.0 {
                    r.warn(&format!(
                        "Docker memory: {gb} GB — 16 GB+ recommended. Some steps may OOM."
                    ));
                } else {
                    r.fail(&format!(
                        "Docker memory: {gb} GB — far too low. Increase to at least 16 GB."
                    ));
                }
            } else {
                r.warn("Could not detect Docker memory allocation");
            }
        } else if have {
            r.fail("Docker daemon is not running. Start Docker Desktop or run: sudo systemctl start docker");
        }

        if command_v("wget") {
            r.pass("wget available (used for downloading reference data)");
        } else if command_v("curl") {
            r.pass("curl available (can substitute for wget in downloads)");
        } else {
            r.fail(
                "Neither wget nor curl found. Install one: apt install wget OR brew install wget",
            );
        }

        let cpus = if is_file("/proc/cpuinfo") {
            match std::fs::read("/proc/cpuinfo") {
                Ok(b) => {
                    let n = lines(&String::from_utf8_lossy(&b))
                        .iter()
                        .filter(|l| l.starts_with("processor"))
                        .count();
                    if n == 0 {
                        "0\n0".to_owned()
                    } else {
                        n.to_string()
                    }
                }
                Err(_) => "0".to_owned(),
            }
        } else if command_v("sysctl") {
            let (ok, o) = capture("sysctl", &["-n", "hw.ncpu"]);
            if ok { o } else { or_zero(o) }
        } else if command_v("nproc") {
            let (ok, o) = capture("nproc", &[]);
            if ok { o } else { or_zero(o) }
        } else {
            "0".to_owned()
        };
        match test_int(&cpus) {
            Some(n) if n > 0 => {
                if n >= 16 {
                    r.pass(&format!("CPU cores: {n} (excellent for parallel steps)"));
                } else if n >= 4 {
                    r.pass(&format!(
                        "CPU cores: {n} (adequate; 16+ recommended for faster runs)"
                    ));
                } else {
                    r.warn(&format!(
                        "CPU cores: {n} — minimum is 4. Pipeline will be very slow."
                    ));
                }
            }
            _ => r.info("Could not detect CPU count"),
        }

        let ram_kb = if is_file("/proc/meminfo") {
            match std::fs::read("/proc/meminfo") {
                Ok(b) => lines(&String::from_utf8_lossy(&b))
                    .iter()
                    .filter(|l| l.starts_with("MemTotal:"))
                    .map(|l| l.split_whitespace().nth(1).unwrap_or("").to_owned())
                    .collect::<Vec<_>>()
                    .join("\n"),
                Err(_) => "0".to_owned(),
            }
        } else if command_v("sysctl") {
            let (ok, o) = capture("sysctl", &["-n", "hw.memsize"]);
            let bytes = if ok { o } else { or_zero(o) };
            (test_int(&bytes).unwrap_or(0) / 1024).to_string()
        } else {
            "0".to_owned()
        };
        match test_int(&ram_kb) {
            Some(kb) if kb > 0 => {
                let gb: i64 = fixed(kb as f64 / 1_048_576.0, 0).parse().unwrap_or(0);
                if gb >= 32 {
                    r.pass(&format!(
                        "System RAM: {gb} GB (can run multiple steps in parallel)"
                    ));
                } else if gb >= 16 {
                    r.pass(&format!(
                        "System RAM: {gb} GB (adequate; 32 GB+ allows more parallelism)"
                    ));
                } else {
                    r.warn(&format!("System RAM: {gb} GB — 16 GB minimum, 32 GB recommended. Reduce --memory flags in scripts."));
                }
            }
            _ => r.info("Could not detect system RAM"),
        }
        // Docker Desktop and colima run the containers in a VM. run-all.sh
        // then sizes the run to the VM, by the same rule (whole GB).
        let engine_gb = if digits(&engine_mem) {
            engine_mem.parse::<u64>().unwrap_or(0) / 1_073_741_824
        } else {
            0
        };
        let host_gb = if digits(&ram_kb) {
            ram_kb.parse::<u64>().unwrap_or(0) / 1_048_576
        } else {
            0
        };
        if engine_gb >= 1 && host_gb >= 1 && engine_gb < host_gb {
            r.info(&format!("Docker has {engine_gb} GB of this host's {host_gb} GB: run-all.sh sizes the tasks to Docker's {engine_gb} GB, not the host's"));
        }

        let free_gb = self.disk(r);

        let arch = self.uname("-m");
        let os = self.uname("-s");
        if arch == "arm64" || arch == "aarch64" {
            r.warn(&format!("Architecture: {arch} ({os}) — most pipeline Docker images are amd64 only. Expect those to run 2-5x slower under emulation."));
        } else {
            r.info(&format!("Architecture: {arch} ({os})"));
        }
        free_gb
    }

    /// `$(uname ARG 2>/dev/null || echo unknown)`.
    fn uname(&self, arg: &str) -> String {
        let (ok, o) = capture("uname", &[arg]);
        if ok {
            o
        } else if o.is_empty() {
            "unknown".to_owned()
        } else {
            format!("{o}\nunknown")
        }
    }

    /// The engine's memory in bytes, as `$(ENGINE info --format
    /// '{{.MemTotal}}' 2>/dev/null || echo 0)` gives it; podman keeps it in
    /// `{{.Host.MemTotal}}`, asked when the first gives no number.
    fn engine_memory(&self) -> String {
        let ask = |f: &str| {
            let (ok, o) = self.engine.capture(&["info", "--format", f]);
            if ok { o } else { or_zero(o) }
        };
        let m = ask("{{.MemTotal}}");
        if test_int(&m).is_some_and(|n| n > 0) {
            return m;
        }
        let h = ask("{{.Host.MemTotal}}");
        if test_int(&h).is_some_and(|n| n > 0) {
            h
        } else {
            m
        }
    }

    /// The work directory of SAMPLE's run: NXF_WORK, which run-all.sh sets
    /// (a relative one under SAMPLE/nextflow), else SAMPLE/nextflow/work.
    fn work_dir(&self) -> String {
        let (g, sm) = (&self.genome, &self.sample);
        match self.vars.nonempty("NXF_WORK") {
            Some(w) if w.starts_with('/') => w.to_owned(),
            Some(w) if !sm.is_empty() && !g.is_empty() => format!("{g}/{sm}/nextflow/{w}"),
            Some(w) => w.to_owned(),
            None if !sm.is_empty() && !g.is_empty() => format!("{g}/{sm}/nextflow/work"),
            None => String::new(),
        }
    }

    /// The work directory when the run is a resume of an earlier one, ""
    /// otherwise: it holds tasks, SAMPLE's launch folder holds a cache of an
    /// earlier session, and its newest log names this work directory.
    fn resume_work(&self, work: &str) -> String {
        let (g, sm) = (&self.genome, &self.sample);
        if sm.is_empty() || g.is_empty() || work.is_empty() || !is_dir(work) {
            return String::new();
        }
        if std::fs::read_dir(work).map_or(true, |mut d| d.next().is_none()) {
            return String::new();
        }
        let cache = format!("{g}/{sm}/nextflow/.nextflow/cache");
        let session = std::fs::read_dir(&cache).ok().is_some_and(|d| {
            d.flatten().any(|e| {
                !e.file_name().to_string_lossy().starts_with('.')
                    && std::fs::symlink_metadata(e.path().join("db")).is_ok()
            })
        });
        if !session {
            return String::new();
        }
        let last = nf_last_work_dir(&format!("{g}/{sm}/nextflow/.nextflow.log"));
        if !last.is_empty() && physical(&last) == physical(work) {
            work.to_owned()
        } else {
            String::new()
        }
    }

    /// The free disk lines: GENOME_DIR's disk (or the working folder's), and
    /// the work directory's when it is on another one.
    fn disk(&self, r: &mut Report) -> Option<i64> {
        let check_dir = if self.genome.is_empty() {
            std::env::current_dir()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_default()
        } else {
            self.genome.clone()
        };
        let work = self.work_dir();
        let resume = self.resume_work(&work);
        let mut told = false;
        let mut free_gb = None;
        if command_v("df") {
            let mut mount = None;
            if let Some((kb, m)) = disk_of(&check_dir)
                && let Some(k) = test_int(&kb).filter(|k| *k >= 0)
            {
                free_gb = fixed(k as f64 / 1_048_576.0, 0).parse::<i64>().ok();
                check_free(r, &kb, &m, "Free disk space", &resume, &mut told);
                mount = Some(m);
            }
            if !work.is_empty()
                && let Some((kb, m)) = disk_of(&work)
                && test_int(&kb).is_some_and(|k| k >= 0)
                && Some(&m) != mount.as_ref()
            {
                check_free(
                    r,
                    &kb,
                    &m,
                    &format!("Free disk space for Nextflow's work directory {work}"),
                    &resume,
                    &mut told,
                );
            }
        }
        free_gb
    }
}

/// `nf_last_work_dir LOG`: the work directory the last
/// "nextflow.Session - Work-dir: PATH [fs]" line of LOG names, "" when none.
pub fn nf_last_work_dir(log: &str) -> String {
    let marker = " nextflow.Session - Work-dir: ";
    let mut last = String::new();
    for l in lines(&read_text(log)) {
        // sed's `^.* nextflow\.Session - Work-dir: \(.*\) \[[^]]*\]$`: the
        // last marker whose rest ends in one bracketed group.
        let mut from = l.len();
        while let Some(i) = l[..from].rfind(marker) {
            let rest = &l[i + marker.len()..];
            if let Some(path) = bracket_tail(rest) {
                last = path.to_owned();
                break;
            }
            if i == 0 {
                break;
            }
            from = i;
        }
    }
    last
}

/// `\(.*\) \[[^]]*\]$` of REST: the text before its last ` [` when the line
/// ends in `]` with no `]` between them.
fn bracket_tail(rest: &str) -> Option<&str> {
    let body = rest.strip_suffix(']')?;
    let i = body.rfind(" [")?;
    (!body[i + 2..].contains(']')).then(|| &body[..i])
}

/// disk_of DIR: the free KB and the mount point of the disk that holds DIR,
/// or its nearest parent that exists (`df -Pk`).
fn disk_of(dir: &str) -> Option<(String, String)> {
    let mut d = dir.to_owned();
    while !exists(&d) && d != "/" && d != "." {
        d = dirname(&d);
    }
    let (_, out) = capture("df", &["-Pk", &d]);
    let l = lines(&out).get(1).map(|l| (*l).to_owned())?;
    let f: Vec<&str> = l.split_whitespace().collect();
    let free = f.get(3).copied().unwrap_or("").to_owned();
    let mount = if f.len() > 5 {
        f[5..].join(" ")
    } else {
        String::new()
    };
    let line = format!("{free} {mount}");
    // read -r FREE_KB CHECK_MOUNT
    let mut it = line.trim().splitn(2, [' ', '\t']);
    let kb = it.next().unwrap_or("").to_owned();
    let m = it.next().unwrap_or("").trim().to_owned();
    if kb.is_empty() { None } else { Some((kb, m)) }
}

/// check_free FREE_KB MOUNT WHAT: the 200 GB minimum and the 500 GB advice.
fn check_free(r: &mut Report, kb: &str, mount: &str, what: &str, resume: &str, told: &mut bool) {
    let k: f64 = kb.trim().parse().unwrap_or(0.0);
    let gb: i64 = fixed(k / 1_048_576.0, 0).parse().unwrap_or(0);
    if gb >= 500 {
        r.pass(&format!("{what}: {gb} GB in {mount}"));
    } else if gb >= 200 {
        r.warn(&format!(
            "{what}: {gb} GB — 500 GB+ recommended for full pipeline per sample"
        ));
    } else if !resume.is_empty() && gb > 0 {
        // A full disk (0 GB) fails even on a resume.
        r.warn(&format!(
            "{what}: {gb} GB, under the 200 GB minimum, but {resume} exists: a resume of an earlier run."
        ));
        if !*told {
            r.say("       -resume reuses every task that finished, so it writes only what the unfinished tasks");
            r.say("       write. A task that must start over needs its full space again: a redone alignment");
            r.say("       writes the BAM (80-120 GB) into work/ and a copy into aligned/. Deleting work/ frees");
            r.say("       its space but makes the next run start from the beginning.");
            *told = true;
        }
    } else {
        r.fail(&format!(
            "{what}: {gb} GB, under the 200 GB minimum (500 GB+ recommended per sample)."
        ));
    }
}

/// The files under DIR named *.bai, *.crai, *.csi or *.tbi, not under
/// SKIP, in the order `find DIR -path SKIP -prune -o -type f ... -print`
/// gives them.
fn walk(dir: &Path, skip: &str, out: &mut Vec<String>) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    for e in rd.flatten() {
        let p = e.path();
        let ps = p.to_string_lossy().into_owned();
        if ps == skip {
            continue;
        }
        let Ok(ft) = e.file_type() else { continue };
        if ft.is_dir() {
            walk(&p, skip, out);
        } else if ft.is_file() {
            let n = e.file_name().to_string_lossy().into_owned();
            if [".bai", ".crai", ".csi", ".tbi"]
                .iter()
                .any(|x| n.ends_with(x))
            {
                out.push(ps);
            }
        }
    }
}

/// A word of an aligner's command line that names reads or alignments
/// (.fastq, .fq, .bam, .cram, .sam, gzip or bzip2 compressed or not).
fn read_file_word(x: &str) -> bool {
    let x = x
        .strip_suffix(".gz")
        .or_else(|| x.strip_suffix(".bz2"))
        .unwrap_or(x);
    [".fastq", ".fq", ".bam", ".cram", ".sam"]
        .iter()
        .any(|e| x.ends_with(e))
}

/// awk's `sub(/.*\//, "", b)`: the text after the last slash.
fn basename_any(x: &str) -> &str {
    x.rsplit('/').next().unwrap_or(x)
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

/// `read -r n kept <<< "$counts"`: the first two words of the first line,
/// the second keeping the rest of the line.
fn vcf_counts(counts: &str) -> (Option<String>, Option<String>) {
    let first = counts.split('\n').next().unwrap_or("");
    let t = first.trim_matches(|c| c == ' ' || c == '\t');
    if t.is_empty() {
        return (None, None);
    }
    let mut it = t.splitn(2, [' ', '\t']);
    let n = it.next().map(str::to_owned);
    let kept = it
        .next()
        .map(|k| k.trim_matches(|c| c == ' ' || c == '\t').to_owned());
    (n, kept)
}

impl Report<'_> {
    /// The verdict of check_vcf_ref_bases from N records read and KEPT
    /// records whose REF matched.
    pub fn vcf_ref_bases_verdict(&mut self, n: Option<&str>, kept: Option<&str>) {
        let (Some(n), Some(kept)) = (n.filter(|x| digits(x)), kept.filter(|x| digits(x))) else {
            self.warn("Could not compare the VCF's REF bases with the reference (bcftools norm failed). Check the build yourself:");
            self.say("       docs/vendor-guide.md, Genome Build.");
            return;
        };
        let n: i64 = n.parse().unwrap_or(0);
        let kept: i64 = kept.parse().unwrap_or(0);
        if n == 0 {
            self.info("VCF has no records to compare with the reference");
            return;
        }
        let bad = n - kept;
        if bad == 0 {
            self.pass(&format!(
                "VCF REF bases match the reference in all of the first {n} records"
            ));
        } else if bad * 20 > n {
            self.fail(&format!(
                "VCF REF bases differ from the reference in {bad} of the first {n} records: the VCF is not on GRCh38"
            ));
            self.say("       The pipeline needs GRCh38, and every position would be read on the wrong build. Get a GRCh38");
            self.say("       VCF from your provider, or start from the FASTQ (docs/vendor-guide.md, Genome Build).");
        } else {
            self.warn(&format!(
                "VCF REF bases differ from the reference in {bad} of the first {n} records"
            ));
        }
    }
}

/// The summary lines and their colours.
fn summary(r: &mut Report, missing_images: usize) {
    r.header("Summary");
    r.say("");
    let c = r.c.clone();
    let (f, w) = (r.failures, r.warnings);
    if f == 0 && w == 0 {
        r.say(&format!(
            "  {}{}All checks passed.{} Your setup is ready to run the pipeline.",
            c.green, c.bold, c.reset
        ));
    } else if f == 0 {
        r.say(&format!(
            "  {}{}All critical checks passed{} with {}{w} warning(s){}.",
            c.green, c.bold, c.reset, c.yellow, c.reset
        ));
        r.say("  The pipeline can run, but review the warnings above for best results.");
    } else {
        r.say(&format!(
            "  {}{}{f} critical issue(s){} and {}{w} warning(s){} found.",
            c.red, c.bold, c.reset, c.yellow, c.reset
        ));
        r.say(&format!(
            "  Fix the {}[FAIL]{} items above before running the pipeline.",
            c.red, c.reset
        ));
    }
    if missing_images > 0 {
        r.say("");
        r.say("  To pull all missing Docker images (~10-15 GB total):");
        r.say("    ./scripts/setup.sh --pull-only");
    }
    r.say("");
    r.say("  Documentation:  docs/00-reference-setup.md   (download reference data)");
    r.say("                  docs/hardware-requirements.md (detailed requirements)");
    r.say("                  docs/vendor-guide.md          (data format help)");
    r.say("");
}

impl Check {
    // --- 2. Environment Variables -----------------------------------------------------------

    fn environment(&self, r: &mut Report) {
        r.header("Environment Variables");
        let g = &self.genome;
        if g.is_empty() {
            r.fail("GENOME_DIR is not set. Export it before running the pipeline:");
            r.say("       export GENOME_DIR=/path/to/your/data");
            return;
        }
        r.pass(&format!("GENOME_DIR is set: {g}"));
        if is_dir(g) {
            r.pass("GENOME_DIR directory exists");
        } else {
            r.fail(&format!("GENOME_DIR directory does not exist: {g}"));
            r.say(&format!("       Create it with: mkdir -p \"{g}\""));
        }
    }

    // --- 3. Reference Data ------------------------------------------------------------------

    fn reference(&self, r: &mut Report) {
        r.header("Reference Data");
        let g = &self.genome;
        if g.is_empty() {
            r.info("Skipping reference data checks (GENOME_DIR not set)");
            return;
        }
        // Every installed path a check below reads, asked of the table first,
        // whichever branch then reads it.
        let qc_rows = components::by_key("sample-qc-data").map_or(&[][..], |c| c.installed);
        let qc: Vec<(String, &str)> = qc_rows
            .iter()
            .map(|t| (self.table("sample-qc-data", t), *t))
            .collect();
        let vep_info = self.table(
            "vep-cache",
            "vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt",
        );
        let clinical_vep = self.table(
            "clinical-data",
            "vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt",
        );
        debug_assert_eq!(vep_info, clinical_vep);
        let pcgr_dir = self.table("clinical-data", "pcgr_data/${PCGR_DATA_BUNDLE}/data");
        let annotsv_genes = self.table(
            "annotsv-data",
            "annotsv_annotations/Annotations_Human/Genes/GRCh38",
        );
        let mane_tx = self.table("annotsv-data", "annotsv_annotations/mane_select_refseq.txt");
        let mane_release = self.table(
            "annotsv-data",
            "annotsv_annotations/mane_select_refseq.release",
        );
        let cyrius = self.table("cyrius", "tools/cyrius-${CYRIUS_VERSION}/INSTALLED");
        let parascopy = self.table(
            "parascopy-data",
            "reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz",
        );
        let kir = self.table("kir-data", "kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat");
        let consumer = self.table(
            "consumer-sites",
            "reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz",
        );
        let panel = self.table(
            "ancestry-panel",
            "reference/pgsc_calc/${PGSC_PANEL}.tar.zst",
        );
        let fine_done = self.table(
            "ancestry-fine-panel",
            "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/model/DONE",
        );

        // --- GRCh38 FASTA ---
        let fasta = &self.ref_fasta;
        if is_file(fasta) {
            let size = std::fs::metadata(fasta).map_or(0, |m| m.len());
            let gb = fixed(size as f64 / 1_073_741_824.0, 1);
            if size > 3_000_000_000 {
                r.pass(&format!("GRCh38 FASTA: {gb} GB"));
            } else {
                r.fail(&format!("GRCh38 FASTA exists but is too small ({gb} GB). Expected > 3 GB. Re-download it."));
                r.say("       See docs/00-reference-setup.md for download instructions.");
            }
        } else {
            r.fail(&format!("GRCh38 FASTA not found at: {fasta}"));
            r.say(&format!("       Download it: ./scripts/setup.sh {g}"));
        }

        // --- FASTA index (.fai) ---
        let fai = format!("{fasta}.fai");
        if is_file(&fai) {
            r.pass("FASTA index (.fai) present");
        } else {
            r.fail(&format!("FASTA index not found at: {fai}"));
            r.say(&format!("       Download it: ./scripts/setup.sh {g}"));
        }

        // --- No ALT or HLA contigs ---
        if is_file(&fai) {
            let text = read_text(&fai);
            let alt: Vec<&str> = lines(&text)
                .iter()
                .map(|l| l.split('\t').next().unwrap_or(""))
                .filter(|n| n.ends_with("_alt") || n.starts_with("HLA-"))
                .collect();
            if alt.is_empty() {
                r.pass(&format!(
                    "Reference has no ALT or HLA contigs ({} sequences)",
                    count_nonempty(&text)
                ));
            } else if self.env("ALLOW_ALT_REFERENCE") == "true" {
                r.warn(&format!(
                    "Reference has {} ALT/HLA contigs (allowed by ALLOW_ALT_REFERENCE=true).",
                    alt.len()
                ));
                r.say("       Reads that match a primary locus and its ALT copy equally get MAPQ 0, so depth");
                r.say("       and calls thin at CYP2D6, the MHC and KIR. See docs/realignment.md.");
            } else {
                r.fail(&format!(
                    "Reference has {} ALT/HLA contigs (first: {}).",
                    alt.len(),
                    alt[0]
                ));
                r.say("       No step aligns ALT-aware, so reads that match a primary locus and its ALT copy");
                r.say("       equally get MAPQ 0 and callers ignore them: depth and calls thin at CYP2D6, the MHC");
                r.say(&format!(
                    "       and KIR. Use the default no-ALT reference (./scripts/setup.sh {g}),"
                ));
                r.say("       or set ALLOW_ALT_REFERENCE=true to keep this one on purpose. See docs/realignment.md.");
            }
        }

        self.clinvar(r);

        // --- Small pinned data files (setup.sh installs them) ---
        let hla = self.v("HLA_DB_RELEASE");
        let gencode = self.v("GENCODE_RELEASE");
        let asset = self.root.join("assets/reference/cytoBand.hg38.txt");
        for name in DATA_FILES {
            let (what, path) = match name {
                "delly_exclude" => (
                    "Delly exclude map (step 19 runs without it, slower and with more artefacts)".to_owned(),
                    format!("{g}/reference/delly_human.hg38.excl.tsv"),
                ),
                "cytoband" => (
                    "GRCh38 chromosome bands (step 10 installs them from assets/reference when missing)".to_owned(),
                    format!("{g}/reference/cytoBand.hg38.txt"),
                ),
                "hla_dat" => (
                    format!("IPD-IMGT/HLA {hla} (step 08 is skipped without it)"),
                    format!("{g}/hla/IPD-IMGT-HLA_{hla}/hla.dat"),
                ),
                "gencode_genes" => (
                    format!("GENCODE {gencode} gene coordinates (step 08 is skipped without them)"),
                    format!("{g}/reference/gencode.v{gencode}.basic.genes.gtf"),
                ),
                _ => (
                    format!("GENCODE {gencode} exons (steps 06 and 17 skip their cDNA read check without them)"),
                    format!("{g}/reference/gencode.v{gencode}.basic.exons.bed"),
                ),
            };
            let short = what.split(" (").next().unwrap_or(&what).to_owned();
            if s(&path) {
                let same = match (std::fs::read(&path), std::fs::read(&asset)) {
                    (Ok(a), Ok(b)) => a == b,
                    _ => false,
                };
                if name == "cytoband" && !same {
                    r.fail(&format!(
                        "{short} (step 10): {path} is not the copy in assets/reference (sha256 {}, expected {})",
                        sha256(Path::new(&path)),
                        sha256(&asset)
                    ));
                    r.say(&format!(
                        "       Delete it and reinstall: ./scripts/setup.sh {g}"
                    ));
                } else {
                    r.pass(&format!("{short}: present"));
                }
            } else {
                r.warn(&format!("{what} not found at: {path}"));
                r.say(&format!("       Install it: ./scripts/setup.sh {g}"));
            }
        }

        // --- gnomAD v3.1's chrM sites, for step 20's annotation ---
        let url = self.v("GNOMAD_MITO_URL");
        let mito = format!("{g}/gnomad_mito/{}", url.rsplit('/').next().unwrap_or(&url));
        if !s(&mito) {
            r.warn(&format!("gnomAD v3.1 mitochondrial sites (step 20 leaves its gnomAD columns empty) not found at: {mito}"));
            r.say(&format!("       Install it: ./scripts/setup.sh {g}"));
        } else if sha256(Path::new(&mito)) != self.v("GNOMAD_MITO_SHA256") {
            r.fail(&format!("gnomAD v3.1 mitochondrial sites (step 20): {mito} is not the pinned file (sha256 differs from GNOMAD_MITO_SHA256)"));
            r.say(&format!(
                "       Delete it and reinstall: ./scripts/setup.sh {g}"
            ));
        } else {
            r.pass("gnomAD v3.1 mitochondrial sites (step 20): present");
        }

        // --- somalier's sites and VerifyBamID2's panel, for step 33 ---
        let missing: String = qc
            .iter()
            .filter(|(p, _)| !s(p))
            .map(|(_, t)| format!(" {t}"))
            .collect();
        if missing.is_empty() {
            r.pass(&format!(
                "somalier sites and VerifyBamID2 panel (step {}): present",
                Self::steps("sample-qc-data")
            ));
        } else {
            r.warn(&format!(
                "Step {} (sample identity and contamination) data not found:{missing}",
                Self::steps("sample-qc-data")
            ));
            r.say(&format!(
                "       Install it: ./scripts/setup.sh {} {g}",
                Self::flag("sample-qc-data")
            ));
        }

        // --- Clinical data installed by default by setup.sh ---
        let vep_rel = self.v("VEP_CACHE_RELEASE");
        let pcgr_rel = self.v("PCGR_VEP_CACHE_RELEASE");
        let vep_release_dir = dirname(&vep_info);
        let vep_dir = dirname(&vep_release_dir);
        if is_file(&vep_info) {
            r.pass(&format!("VEP {vep_rel} cache (step 13): present"));
        } else {
            r.warn(&format!(
                "VEP {vep_rel} cache not found at: {vep_release_dir}"
            ));
            r.say("       Missing cache skips steps 13 (VEP), 30 (vcfanno), 23 (clinical filter) and 31 (slivar).");
            r.say(&format!(
                "       Install it: ./scripts/setup.sh --clinical-data \"{g}\""
            ));
            r.say("       See docs/00-reference-setup.md for full instructions.");
        }
        if is_file(&format!("{vep_dir}/{pcgr_rel}_GRCh38/info.txt")) {
            r.pass(&format!("VEP {pcgr_rel} cache (step 17 CPSR): present"));
        } else {
            r.warn(&format!(
                "VEP {pcgr_rel} cache not found at: {vep_dir}/{pcgr_rel}_GRCh38"
            ));
            r.say(&format!(
                "       Required for step 17 (CPSR). {} needs VEP {pcgr_rel}, separate from step 13's VEP {vep_rel}.",
                self.v("PCGR_IMAGE")
            ));
            r.say(&format!("       Missing cache skips step 17. Install it: ./scripts/setup.sh --clinical-data \"{g}\""));
            r.say("       See docs/17-cpsr.md for full instructions.");
        }

        // --- PCGR data bundle (optional) ---
        if is_dir(&pcgr_dir) {
            r.pass("PCGR/CPSR ref data bundle: present");
        } else {
            r.warn(&format!("PCGR ref data bundle not found at: {pcgr_dir}"));
            r.say("       Missing bundle skips step 17 (CPSR cancer predisposition).");
            r.say(&format!(
                "       Install it: ./scripts/setup.sh --clinical-data \"{g}\""
            ));
            r.say("       See docs/17-cpsr.md for full instructions.");
        }

        // --- AnnotSV annotation data (step 5) ---
        let annotsv_dir = dirname(&mane_tx);
        if is_dir(&annotsv_genes) {
            r.pass("AnnotSV annotation data (step 5): present");
            let sv = self.v("GNOMAD_SV_VERSION");
            let stamp = format!(
                "{annotsv_dir}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed"
            );
            if cat(&stamp) == format!("gnomAD-SV {sv} sha256 {}", self.v("GNOMAD_SV_SHA256")) {
                r.pass(&format!(
                    "AnnotSV benign set (step 5): gnomAD-SV {sv} common SVs added"
                ));
            } else {
                r.fail(&format!("AnnotSV benign set (step 5): gnomAD-SV {sv} not added ({stamp} missing or of another release)"));
                r.say(&format!("       Without it AnnotSV grades common deletions as if no one carried them. Add it with: ./scripts/setup.sh --annotsv-data {g}"));
            }
            let mane = self.v("MANE_RELEASE");
            if s(&mane_tx)
                && cat(&mane_release)
                    == format!("MANE {mane} sha256 {}", self.v("MANE_SUMMARY_SHA256"))
            {
                let n = std::fs::read(&mane_tx)
                    .map(|b| b.iter().filter(|c| **c == b'\n').count())
                    .unwrap_or(0);
                r.pass(&format!(
                    "AnnotSV transcripts (step 5): MANE {mane} Select list, {n} transcripts"
                ));
            } else {
                r.fail(&format!("AnnotSV transcripts (step 5): MANE {mane} Select list missing or of another release ({mane_tx})"));
                r.say(&format!("       Without it the Tx, Frameshift and exon columns follow whichever transcript overlaps most. Add it with: ./scripts/setup.sh --annotsv-data {g}"));
            }
        } else {
            r.warn(&format!(
                "AnnotSV annotation data not found at: {annotsv_genes}"
            ));
            r.say(&format!("       Step 5 (AnnotSV) will be skipped. Download it (~5.3 GB, ~20 GB unpacked) with: ./scripts/setup.sh {g}"));
        }

        self.annotation_scores(r);

        // --- Opt-in steps: Cyrius (21), Parascopy (35), IPD-KIR (08 with KIR=true), the raw-data export (40) ---
        let cv = self.v("CYRIUS_VERSION");
        let lock = sha256(&self.root.join("scripts/cyrius-constraints.txt"));
        let py = self.v(Self::image_var("cyrius", "PYTHON_IMAGE"));
        if cat(&cyrius) == format!("python={py} lock={lock}") {
            r.pass(&format!(
                "Cyrius {cv} (opt-in step {}): installed",
                Self::steps("cyrius")
            ));
        } else {
            r.info(&format!(
                "Cyrius {cv} not installed (opt-in step {}; non-commercial licence): ./scripts/setup.sh {} {g}",
                Self::steps("cyrius"),
                Self::flag("cyrius")
            ));
        }
        let pv = self.v("PARASCOPY_DATA_VERSION");
        if s(&parascopy) {
            r.pass(&format!(
                "Parascopy {pv} homology table and models (opt-in step {}): present",
                Self::steps("parascopy-data")
            ));
        } else {
            r.info(&format!(
                "Parascopy data not installed (opt-in step {}): ./scripts/setup.sh {} {g}",
                Self::steps("parascopy-data"),
                Self::flag("parascopy-data")
            ));
        }
        let kv = self.v("KIR_DB_RELEASE");
        if s(&kir) {
            r.pass(&format!(
                "IPD-KIR {kv} (step {} with KIR=true): present",
                Self::steps("kir-data")
            ));
        } else {
            r.info(&format!(
                "IPD-KIR {kv} not installed (KIR=true in step {}): ./scripts/setup.sh {} {g}",
                Self::steps("kir-data"),
                Self::flag("kir-data")
            ));
        }
        let build = self.v("CONSUMER_DBSNP_BUILD");
        let flag = Self::flag("consumer-sites");
        let step = Self::steps("consumer-sites");
        if s(&consumer) && gzip_t(&consumer) {
            if consumer_sites_digest(&consumer) == self.v("CONSUMER_SITES_TABLE_SHA256") {
                r.pass(&format!("Site table of the raw-data export (opt-in step {step}, dbSNP {build}): present, the same as the pinned build"));
            } else {
                r.warn(&format!("Site table {consumer} differs from the pinned build (CONSUMER_SITES_TABLE_SHA256): it was built by another version of the pipeline. Step {step} uses it; to rebuild it, remove it and run ./scripts/setup.sh {flag} {g}"));
            }
        } else if exists(&consumer) {
            r.fail(&format!(
                "Site table {consumer} is damaged: remove it and run ./scripts/setup.sh {flag} {g}"
            ));
        } else {
            r.info(&format!("Site table of the raw-data export not built (opt-in step {step}): ./scripts/setup.sh {flag} {g}"));
        }

        // --- Step 25 (pgsc_calc, its scores) and its ancestry panel (step 26) ---
        let (mut n_pgs, mut n_acad) = (0, 0);
        let (mut missing_pgs, mut missing_acad) = (String::new(), String::new());
        for (id, lic) in pgs_rows(&self.root.join("assets/pgs_scores.tsv")) {
            let there = s(&format!("{g}/prs_scores/{id}.txt.gz"));
            if lic == "academic-research" {
                n_acad += 1;
                if !there {
                    missing_acad.push_str(&format!(" {id}"));
                }
                continue;
            }
            n_pgs += 1;
            if !there {
                missing_pgs.push_str(&format!(" {id}"));
            }
        }
        if missing_pgs.is_empty() {
            r.pass(&format!(
                "PGS Catalog scores (step 25): all {n_pgs} of assets/pgs_scores.tsv present"
            ));
        } else {
            r.warn(&format!("PGS Catalog scores missing (step 25 downloads them, run-all.sh skips PRS until then):{missing_pgs}"));
        }
        if n_acad > 0 && missing_acad.is_empty() {
            r.pass(&format!("Academic-only PGS Catalog scores (step 25 with PRS_ACADEMIC=true): all {n_acad} present"));
        } else if n_acad > 0 {
            r.info(&format!(
                "Academic-only PGS Catalog scores not installed (opt-in, PRS_ACADEMIC=true):{missing_acad}: ./scripts/setup.sh {} {g}",
                Self::flag("prs-academic")
            ));
        }
        let pc = self.v("PGSC_CALC_VERSION");
        if is_file(&format!("{g}/tools/pgsc_calc-{pc}/main.nf")) {
            r.pass(&format!("pgsc_calc {pc} (step 25): present"));
        } else {
            r.warn(&format!(
                "pgsc_calc {pc} not installed (step 25 fetches it): ./scripts/setup.sh {g}"
            ));
        }
        let pflag = Self::flag("ancestry-panel");
        let sites = format!(
            "{}_GRCh38_sites.tsv",
            panel.strip_suffix(".tar.zst").unwrap_or(&panel)
        );
        if s(&panel) && s(&sites) {
            r.pass(&format!(
                "Ancestry panel {} (step 26, percentiles in step 25): present",
                self.v("PGSC_PANEL")
            ));
        } else if s(&panel) {
            r.fail(&format!(
                "Ancestry panel {panel} has no site list beside it: ./scripts/setup.sh {pflag} {g}"
            ));
        } else {
            r.info(&format!("Ancestry panel not installed (step 26; step 25 reports raw scores only): ./scripts/setup.sh {pflag} {g}"));
        }
        // Step 26b's panel: model/DONE is the last file setup.sh --ancestry-fine-panel writes.
        let af_panel = dirname(&dirname(&fine_done));
        let fflag = Self::flag("ancestry-fine-panel");
        let fstep = Self::steps("ancestry-fine-panel");
        if s(&fine_done) {
            r.pass(&format!(
                "Fine-grained ancestry panel {} (step {fstep}): present",
                self.v("ANCESTRY_FINE_PANEL")
            ));
        } else if is_dir(&af_panel) {
            r.warn(&format!("Fine-grained ancestry panel {af_panel} is not finished (no model/DONE): ./scripts/setup.sh {fflag} {g}"));
        } else {
            r.info(&format!("Fine-grained ancestry panel not installed (step {fstep}, opt-in): ./scripts/setup.sh {fflag} {g}"));
        }

        // --- GATK sequence dictionary (optional) ---
        let dict = &self.ref_dict;
        if is_file(dict) {
            r.pass("GATK sequence dictionary (.dict): present");
        } else {
            r.warn(&format!("GATK sequence dictionary not found at: {dict}"));
            r.say("       Some tools (GATK Mutect2/step 20) require it. Generate with:");
            r.say(&format!(
                "       docker run --rm -v \"{g}:/genome\" {} \\",
                self.v("GATK_IMAGE")
            ));
            r.say(&format!(
                "         gatk CreateSequenceDictionary -R {}",
                self.ref_fasta_c
            ));
        }

        self.check_reference_content(r);
    }

    /// The ClinVar lines: the two files step 06 needs (a warning when TOOLS
    /// names no clinvar), the conflicting subset, and the release's age.
    fn clinvar(&self, r: &mut Report) {
        let g = &self.genome;
        let clinvar = format!("{g}/clinvar/clinvar_chr.vcf.gz");
        let path = format!("{g}/clinvar/clinvar_pathogenic_chr.vcf.gz");
        let tools = self.env("TOOLS");
        let listed = format!(",{},", tools.replace(' ', "")).contains(",clinvar,");
        if !tools.is_empty() && !listed && (!is_file(&clinvar) || !is_file(&path)) {
            r.warn(&format!("ClinVar not installed (clinvar/clinvar_chr.vcf.gz or clinvar/clinvar_pathogenic_chr.vcf.gz): TOOLS ({tools}) has no clinvar, so step 06 does not run; steps 13 and 20 annotate without ClinVar."));
            r.say(&format!("       To add it: ./scripts/setup.sh {g}"));
        } else {
            if is_file(&clinvar) {
                r.pass("ClinVar VCF (chr-prefixed): present");
            } else if is_file(&format!("{g}/clinvar/clinvar.vcf.gz")) {
                r.warn("ClinVar raw VCF found but chr-prefixed version missing.");
                r.say(&format!(
                    "       Run: ./scripts/setup.sh {g}  OR  see docs/00-reference-setup.md"
                ));
            } else {
                r.fail(&format!("ClinVar VCF not found at: {clinvar}"));
                r.say("       Download and prepare it — see docs/00-reference-setup.md");
            }
            if is_file(&path) {
                r.pass("ClinVar pathogenic subset: present");
            } else if is_file(&clinvar) {
                r.fail(
                    "ClinVar chr-prefixed found but pathogenic subset missing (step 6 will fail).",
                );
                r.say(&format!(
                    "       Run: ./scripts/setup.sh {g}  OR  see docs/00-reference-setup.md"
                ));
            } else {
                r.fail("ClinVar pathogenic subset not found (step 6 will fail).");
                r.say("       Download and prepare it — see docs/00-reference-setup.md");
            }
            let conf = format!("{g}/clinvar/clinvar_conflicting_plp_chr.vcf.gz");
            let conf_tbi = format!("{conf}.tbi");
            if is_file(&conf)
                && is_file(&conf_tbi)
                && !(is_file(&clinvar) && newer(&clinvar, &conf))
            {
                r.pass("ClinVar conflicting subset: present");
            } else if is_file(&path) {
                if is_file(&conf) && is_file(&conf_tbi) {
                    r.info("ClinVar conflicting subset is older than clinvar/clinvar_chr.vcf.gz: step 6 lists conflicting records of an older release.");
                } else {
                    r.info("ClinVar conflicting subset not built: step 6's table of conflicting records some labs call pathogenic says \"not checked\".");
                }
                r.say(&format!("       Build it: ./scripts/setup.sh {g}"));
            }
        }

        // NCBI publishes ClinVar monthly; a copy over 35 days old misses the last release.
        if is_file(&format!("{g}/clinvar/clinvar.vcf.gz")) {
            let rel_file = format!("{g}/clinvar/RELEASE");
            let release = std::fs::read(&rel_file)
                .map(|b| {
                    let t = String::from_utf8_lossy(&b).into_owned();
                    t.split('\n').next().unwrap_or("").to_owned()
                })
                .unwrap_or_default();
            if is_date(&release) {
                match days_since(&release) {
                    Some(age) if age > 35 => {
                        r.warn(&format!(
                            "ClinVar release {release} is {age} days old (over 35)."
                        ));
                        r.say(&format!("       Replace it with the current release: ./scripts/setup.sh --refresh clinvar {g}"));
                    }
                    Some(age) => r.pass(&format!("ClinVar release {release} ({age} days old)")),
                    None => r.warn(&format!(
                        "ClinVar release {release}: could not work out its age on this system"
                    )),
                }
            } else {
                r.warn(&format!(
                    "ClinVar release date unknown (no {g}/clinvar/RELEASE)."
                ));
                r.say(&format!("       Record it with the current release: ./scripts/setup.sh --refresh clinvar {g}"));
            }
        }
    }

    /// The seven score tables of steps 30 and 31, in annotations/ (where
    /// setup.sh --annotation-scores writes them).
    fn annotation_scores(&self, r: &mut Report) {
        let d = format!("{}/annotations", self.genome);
        let both = |f: &str| is_file(&format!("{d}/{f}")) && is_file(&format!("{d}/{f}.tbi"));
        let total = 7;
        let mut count = 0;
        let mut missing: Vec<String> = Vec::new();
        let mut one = |ok: bool, what: String| {
            if ok {
                count += 1;
            } else {
                missing.push(what);
            }
        };
        one(
            both("whole_genome_SNVs.tsv.gz"),
            "CADD SNVs (whole_genome_SNVs.tsv.gz + .tbi)".to_owned(),
        );
        one(
            both("gnomad.genomes.r4.0.indel.tsv.gz"),
            "CADD indels (gnomad.genomes.r4.0.indel.tsv.gz + .tbi)".to_owned(),
        );
        for kind in ["snv", "indel"] {
            let ok = both(&format!("spliceai_scores.raw.{kind}.hg38.vcf.gz"))
                || both(&format!("spliceai_scores.masked.{kind}.hg38.vcf.gz"));
            one(
                ok,
                format!(
                    "SpliceAI {kind}s (spliceai_scores.raw.{kind}.hg38.vcf.gz or spliceai_scores.masked.{kind}.hg38.vcf.gz, + .tbi)"
                ),
            );
        }
        one(
            both("revel_grch38.tsv.gz"),
            "REVEL (revel_grch38.tsv.gz + .tbi)".to_owned(),
        );
        one(
            both("AlphaMissense_hg38.tsv.gz"),
            "AlphaMissense (AlphaMissense_hg38.tsv.gz + .tbi)".to_owned(),
        );
        one(
            is_file(&format!("{d}/gnomad_v4.1_constraint.tsv")),
            "gnomAD constraint (gnomad_v4.1_constraint.tsv)".to_owned(),
        );
        if count == total {
            r.pass(&format!("Annotation databases: all {total} present (CADD SNVs+indels, SpliceAI SNVs+indels, REVEL, AlphaMissense, gnomAD constraint)"));
        } else if count > 0 {
            r.warn(&format!("Annotation databases: {count}/{total} present. Steps 30-31 will degrade gracefully for missing tracks."));
            for m in &missing {
                r.say(&format!("         Missing: {m}"));
            }
            r.say("       See docs/00-reference-setup.md for download instructions.");
        } else {
            r.info("Annotation databases not downloaded (steps 30-31 will be skipped). See docs/00-reference-setup.md");
        }
    }

    // --- 4. Docker Images -------------------------------------------------------------------

    /// The image lines; the number of images not pulled.
    fn images(&self, r: &mut Report) -> usize {
        r.header("Docker Images");
        // The opt-in images and the data that makes them matter, from their rows.
        let opts: Vec<(&str, &str, String, &str)> = [
            (
                "cyrius",
                "PYTHON_IMAGE",
                "tools/cyrius-${CYRIUS_VERSION}/INSTALLED",
            ),
            (
                "parascopy-data",
                "PARASCOPY_IMAGE",
                "reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz",
            ),
            (
                "yleaf-data",
                "YLEAF_IMAGE",
                "reference/yleaf-${YLEAF_DATA_VERSION}/data/hg38/new_positions.txt",
            ),
            (
                "ancestry-fine-panel",
                "BEAGLE_IMAGE",
                "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/model/DONE",
            ),
        ]
        .iter()
        .map(|(k, var, t)| {
            let data = if self.genome.is_empty() {
                String::new()
            } else {
                self.table(k, t)
            };
            (*k, Self::image_var(k, var), data, Self::flag(k))
        })
        .collect();
        let e = &self.engine;
        if !command_v(&e.cmd) || !e.quiet(&["info"]) {
            r.info("Skipping Docker image checks (Docker not available)");
            return 0;
        }
        let imgs = images::pipeline_images(&self.root).unwrap_or_default();
        let (mut pulled, mut missing) = (0, Vec::new());
        let (mut untagged, mut no_pinned, mut unknown) = (Vec::new(), Vec::new(), Vec::new());
        for img in &imgs {
            if e.has(img) {
                pulled += 1;
                if !images::pinned_tag_present(e, img) {
                    match images::image_is_dangling(e, img) {
                        0 => untagged.push(img.clone()),
                        1 => no_pinned.push(img.clone()),
                        _ => unknown.push(img.clone()),
                    }
                }
            } else {
                missing.push(img.clone());
            }
        }
        if missing.is_empty() {
            r.pass(&format!("All {} Docker images are pulled", imgs.len()));
        } else {
            r.pass(&format!(
                "{pulled}/{} Docker images already pulled",
                imgs.len()
            ));
            r.fail(&format!(
                "{} Docker image(s) not yet pulled:",
                missing.len()
            ));
            for img in &missing {
                match images::pinned_tag(img) {
                    Some(t) => r.say(&format!(
                        "         docker pull {img} && docker tag {img} {t}"
                    )),
                    None => r.say(&format!("         docker pull {img}")),
                }
            }
            r.say("");
            r.say("       Pull all missing images at once (and tag the digest-pinned ones):");
            r.say("         ./scripts/setup.sh --pull-only");
        }
        if !untagged.is_empty() {
            r.warn(&format!("{} digest-pinned image(s) have no local tag, so a dangling-image cleanup (docker image prune) would delete them:", untagged.len()));
            for i in &untagged {
                r.say(&format!("         {i}"));
            }
            r.say("       Tag them with: ./scripts/setup.sh --pull-only");
        }
        if !no_pinned.is_empty() {
            r.info(&format!("{} digest-pinned image(s) lack their localhost/pgp-pinned/ tag; Docker does not list them as dangling, so a dangling-image cleanup keeps them:", no_pinned.len()));
            for i in &no_pinned {
                r.say(&format!("         {i}"));
            }
            r.say("       Add it with: ./scripts/setup.sh --pull-only");
        }
        if !unknown.is_empty() {
            r.warn(&format!("{} digest-pinned image(s) are missing their localhost/pgp-pinned/ tag and Docker did not list its dangling images, so a dangling-image cleanup may delete them:", unknown.len()));
            for i in &unknown {
                r.say(&format!("         {i}"));
            }
            r.say("       Tag them with: ./scripts/setup.sh --pull-only");
        }
        // An opt-in step's image matters once its data or the tool is installed.
        if !self.genome.is_empty() {
            for (k, var, data, flag) in &opts {
                if !s(data) {
                    continue;
                }
                let img = self.v(var);
                let step = Self::steps(k);
                if e.has(&img) {
                    r.pass(&format!(
                        "Opt-in step {step}: its data and image {img} are present"
                    ));
                } else {
                    r.warn(&format!("Opt-in step {step}: its data is installed but its image is not pulled; the step would pull it mid-run"));
                    r.say(&format!(
                        "         ./scripts/setup.sh {flag} {}   (or: docker pull {img})",
                        self.genome
                    ));
                }
            }
        }
        missing.len()
    }

    // --- OwnGenome's og binary --------------------------------------------------------------

    fn og_binary(&self, r: &mut Report) {
        r.header("OwnGenome binary");
        // A single build is enough on x86_64. Check both when both were
        // installed, or on Arm where tools may run in amd64 images.
        let set = self.vars.nonempty("OG_ARCHES");
        let mut arches = set.unwrap_or("x86_64 aarch64").to_owned();
        if set.is_none()
            && self.uname_plain() == "x86_64"
            && !exists(&self.root.join("bin/.og/og-linux-aarch64").to_string_lossy())
        {
            arches = "x86_64".to_owned();
        }
        if !is_executable(&self.root.join("bin/og")) {
            r.fail("bin/og is missing or not executable. Restore it from the checkout and run scripts/setup.sh --og.");
            return;
        }
        let og_sh = self.root.join("rust/og.sh");
        let o = Command::new("bash")
            .args([
                "-c",
                "exec bash \"$0\" check --required 2>&1",
                &og_sh.to_string_lossy(),
            ])
            .env("OG_ARCHES", &arches)
            .stdin(Stdio::null())
            .stderr(Stdio::null())
            .output();
        let (ok, msg) = match o {
            Ok(o) => (
                o.status.success(),
                String::from_utf8_lossy(&o.stdout)
                    .trim_end_matches('\n')
                    .to_owned(),
            ),
            Err(_) => (false, String::new()),
        };
        if ok {
            r.pass(msg.strip_prefix("OK: ").unwrap_or(&msg));
        } else {
            r.fail(msg.strip_prefix("ERROR: ").unwrap_or(&msg));
        }
    }

    /// `$(uname -m)`.
    fn uname_plain(&self) -> String {
        capture("uname", &["-m"]).1
    }

    // --- 5. Sample Data ---------------------------------------------------------------------

    fn sample_data(&self, r: &mut Report, free_gb: Option<i64>) {
        let (g, sm) = (&self.genome, &self.sample);
        r.header(&format!("Sample Data: {sm}"));
        if g.is_empty() {
            r.info("Skipping sample checks (GENOME_DIR not set)");
            return;
        }
        if !is_dir(g) {
            r.info("Skipping sample checks (GENOME_DIR does not exist)");
            return;
        }
        let sd = format!("{g}/{sm}");
        let fq = format!("{sd}/fastq");
        let (mut has_ora, mut has_fastq, mut has_bam, mut has_vcf) = (false, false, false, false);

        let ora = if is_dir(&fq) {
            count_in(&fq, |n| n.ends_with(".ora"))
        } else {
            0
        };
        if ora > 0 {
            has_ora = true;
            r.pass(&format!("ORA files found: {ora} file(s) in {fq}/"));
        }

        let r1 = format!("{fq}/{sm}_R1.fastq.gz");
        let r2 = format!("{fq}/{sm}_R2.fastq.gz");
        if is_file(&r1) && is_file(&r2) {
            has_fastq = true;
            let size = std::fs::metadata(&r1).map_or(0, |m| m.len());
            r.pass(&format!(
                "FASTQ files found: R1 ({} GB) + R2",
                fixed(size as f64 / 1_073_741_824.0, 1)
            ));
        } else if is_file(&r1) {
            r.warn(&format!(
                "Only R1 found at {r1} — missing R2. Pipeline expects paired-end reads."
            ));
        } else {
            let n = if is_dir(&fq) {
                count_in(&fq, |n| {
                    [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
                        .iter()
                        .any(|e| n.ends_with(e))
                })
            } else {
                0
            };
            if n > 0 {
                r.warn(&format!("Found {n} FASTQ file(s) in {fq}/ but not named {sm}_R1.fastq.gz / {sm}_R2.fastq.gz"));
                r.say("       Rename or symlink them:");
                r.say(&format!(
                    "         ln -s your_file_R1.fastq.gz {fq}/{sm}_R1.fastq.gz"
                ));
                r.say(&format!(
                    "         ln -s your_file_R2.fastq.gz {fq}/{sm}_R2.fastq.gz"
                ));
            }
        }

        let samtools = self.v("SAMTOOLS_IMAGE");
        let bcftools = self.v("BCFTOOLS_IMAGE");
        let bam = format!("{sd}/aligned/{sm}_sorted.bam");
        let bai = format!("{bam}.bai");
        if is_file(&bam) {
            has_bam = true;
            let size = std::fs::metadata(&bam).map_or(0, |m| m.len());
            r.pass(&format!(
                "BAM file found: {} GB",
                fixed(size as f64 / 1_073_741_824.0, 1)
            ));
            if is_file(&bai) {
                r.pass("BAM index (.bai) present");
            } else {
                r.warn("BAM index not found. Create it before running BAM-dependent steps:");
                r.say(&format!(
                    "       docker run --rm -v \"{g}:/genome\" {samtools} \\"
                ));
                r.say(&format!(
                    "         samtools index /genome/{sm}/aligned/{sm}_sorted.bam"
                ));
            }
        }

        let vcf = format!("{sd}/vcf/{sm}.vcf.gz");
        let tbi = format!("{vcf}.tbi");
        if is_file(&vcf) {
            has_vcf = true;
            let size = std::fs::metadata(&vcf).map_or(0, |m| m.len());
            r.pass(&format!(
                "VCF file found: {} MB",
                fixed(size as f64 / 1_048_576.0, 0)
            ));
            if is_file(&tbi) {
                r.pass("VCF index (.tbi) present");
            } else {
                r.warn("VCF index not found. Create it before running VCF-dependent steps:");
                r.say(&format!(
                    "       docker run --rm -v \"{g}:/genome\" {bcftools} \\"
                ));
                r.say(&format!(
                    "         bcftools index -t /genome/{sm}/vcf/{sm}.vcf.gz"
                ));
            }
        }

        self.check_stale_indexes(r, &sd);

        match self.vars.nonempty("INPUT").unwrap_or("auto") {
            "fastq" => {
                if has_bam {
                    r.info("INPUT=fastq: the existing BAM is not read, so it is not checked");
                }
                if has_vcf {
                    r.info("INPUT=fastq: the existing VCF is not read, so it is not checked");
                }
                has_bam = false;
                has_vcf = false;
            }
            "bam" => {
                if has_vcf {
                    r.info("INPUT=bam: the existing VCF is not read, so it is not checked");
                }
                has_vcf = false;
            }
            "vcf" => {
                if has_bam {
                    r.info("INPUT=vcf: the BAM is not read, so it is not checked");
                }
                has_bam = false;
            }
            _ => {}
        }

        let engine = command_v(&self.engine.cmd);
        let bam_c = self.bam_c();
        if has_bam && engine {
            r.say("");
            r.info("Checking genome build of BAM...");
            let mut len = self.bam_chr1_len(&samtools, &bam_c, |l| l.contains("SN:chr1"));
            if len.is_empty() || len == "0" {
                // Try without chr prefix (hg19 style)
                len = self.bam_chr1_len(&samtools, &bam_c, |l| {
                    l.match_indices("SN:1")
                        .any(|(i, _)| l[i + 4..].chars().next().is_some_and(|c| c.is_whitespace()))
                });
                if !len.is_empty() && len != "0" {
                    r.fail("BAM uses chromosome names WITHOUT 'chr' prefix (hg19/GRCh37 style)");
                    r.say("       This pipeline requires GRCh38 (hg38) with 'chr' prefix.");
                    r.say(
                        "       Extract FASTQ and re-align: samtools fastq -> step 2 (alignment)",
                    );
                    r.say("       See docs/vendor-guide.md for build conversion instructions.");
                }
            } else if len == "248956422" {
                r.pass("BAM genome build: GRCh38 (chr1 length = 248,956,422)");
            } else if len == "249250621" {
                r.fail("BAM genome build: GRCh37/hg19 (chr1 length = 249,250,621)");
                r.say("       This pipeline requires GRCh38. Re-align from FASTQ.");
            } else {
                r.warn(&format!(
                    "BAM chr1 length ({len}) does not match known builds"
                ));
                r.say("       Expected: 248956422 (GRCh38) or 249250621 (GRCh37)");
            }

            // Read group: GATK steps (20, 03a, 29) reject reads without one,
            // and DeepVariant takes the sample name from it.
            let (rc, header) = self.run_in(&[&samtools, "samtools", "view", "-H", &bam_c]);
            if rc == 0 {
                let m5 = self.check_bam_reference(r, &header);
                self.check_bam_masked_reference(r, &header, m5);
                if lines(&header).iter().any(|l| l.starts_with("@RG")) {
                    r.pass("BAM has a read group (@RG)");
                } else {
                    r.warn("BAM header has no @RG read group line. GATK steps (20, 03a, 29) will reject its reads.");
                    r.say("       Add one, then replace the BAM and re-index it:");
                    r.say(&format!(
                        "       docker run --rm -v \"{g}:/genome\" {samtools} \\"
                    ));
                    r.say(&format!("         samtools addreplacerg -r ID:{sm} -r SM:{sm} -r PL:ILLUMINA -r LB:{sm} \\"));
                    r.say(&format!(
                        "         -o /genome/{sm}/aligned/{sm}_sorted.rg.bam /genome/{sm}/aligned/{sm}_sorted.bam"
                    ));
                    r.say(&format!(
                        "       mv \"{sd}/aligned/{sm}_sorted.rg.bam\" \"{bam}\""
                    ));
                    r.say(&format!(
                        "       docker run --rm -v \"{g}:/genome\" {samtools} \\"
                    ));
                    r.say(&format!(
                        "         samtools index /genome/{sm}/aligned/{sm}_sorted.bam"
                    ));
                }
            } else {
                r.fail("Could not read the BAM header (samtools view -H failed): the BAM is unreadable or not a BAM");
            }
            self.check_bam_quickcheck(r);
            self.check_bam_read_count(r, &sd);
        }

        let mut vcf_header = None;
        if has_vcf && engine {
            let vc = format!("/genome/{sm}/vcf/{sm}.vcf.gz");
            let (rc, h) = self.run_in(&[&bcftools, "bcftools", "view", "-h", &vc]);
            if rc == 0 {
                vcf_header = Some(h);
            } else {
                r.fail("Could not read the VCF header (bcftools view -h failed): the VCF is unreadable or not a VCF");
            }
        }

        // Genome build from the header: chr1's ##contig length. When the
        // header cannot tell, compare the REF bases of the first records.
        if has_vcf && let Some(h) = &vcf_header {
            let hl = lines(h);
            let chr1 = hl
                .iter()
                .find(|l| l.starts_with("##contig=<ID=chr1,") || l.starts_with("##contig=<ID=1,"))
                .map(|l| contig_lengths(l))
                .unwrap_or_default();
            if hl.iter().any(|l| l.starts_with("##contig=<ID=1,")) {
                r.fail("VCF uses chromosome names without the 'chr' prefix (1, 2, MT); the pipeline needs chr1, chr2, chrM");
                r.say("       Rename them with bcftools annotate --rename-chrs: docs/vcf-first.md shows the commands.");
            }
            if chr1 == "248956422" {
                r.pass("VCF genome build: GRCh38");
            } else if chr1 == "249250621" {
                r.fail("VCF genome build: GRCh37/hg19 — not compatible with this pipeline");
            } else {
                if !chr1.is_empty() {
                    r.warn(&format!("VCF chr1 length ({chr1}) is neither GRCh38's (248956422) nor GRCh37's (249250621)"));
                } else if hl.iter().any(|l| l.starts_with("##contig=")) {
                    r.warn("VCF header has no chr1 length, so it does not show the genome build");
                } else {
                    r.warn(
                        "VCF header has no ##contig lines, so it does not show the genome build",
                    );
                }
                self.check_vcf_ref_bases(r);
            }

            // One sample per run: a joint-called VCF holds one column per person.
            let samples = vcf_header_samples(h);
            let n = if samples.is_empty() {
                0
            } else {
                samples.split('\n').filter(|x| !x.is_empty()).count()
            };
            if n > 1 {
                let mut names = samples
                    .split('\n')
                    .take(5)
                    .collect::<Vec<_>>()
                    .join(",")
                    .replace(',', ", ");
                if n > 5 {
                    names.push_str(", ...");
                }
                r.fail(&format!(
                    "VCF holds {n} samples ({names}); the pipeline analyses one sample per run"
                ));
                r.say("       Give each sample its own folder with its column only, for example:");
                r.say(&format!("       mkdir -p \"{g}/<name>/vcf\""));
                r.say(&format!(
                    "       docker run --rm -v \"{g}:/genome\" {bcftools} \\"
                ));
                r.say(&format!("         bcftools view -s <name> -a -c 1 -Oz -o /genome/<name>/vcf/<name>.vcf.gz /genome/{sm}/vcf/{sm}.vcf.gz"));
                r.say(&format!(
                    "       docker run --rm -v \"{g}:/genome\" {bcftools} \\"
                ));
                r.say("         bcftools index -t /genome/<name>/vcf/<name>.vcf.gz");
                r.say("       (-a drops the ALT alleles that sample does not carry, -c 1 the sites where it carries none:");
                r.say("       the result is a variant-only VCF)");
            } else if n == 1 {
                r.pass(&format!("VCF holds one sample ({samples})"));
            }
        }

        // A FASTQ start writes far more than a BAM or VCF start.
        if has_fastq
            && !has_bam
            && !has_vcf
            && let Some(free) = free_gb
            && free < 450
        {
            r.warn(&format!("FASTQ start with {free} GB free: alignment alone can write about 450 GB before you delete work/"));
            r.say("       Trimmed reads 60-90 GB, alignment intermediate 119 GB in one observed 30x run, the BAM");
            r.say(&format!(
                "       80-120 GB in {sm}/nextflow/work and its copy in {sm}/aligned. Free more space,"
            ));
            r.say("       or put GENOME_DIR on a larger disk.");
        }

        // Suggest pipeline entry path.
        r.say("");
        let hint = format!("       Whole pipeline:  ./scripts/run-all.sh {sm} <male|female>");
        if has_ora && !has_fastq && !has_bam && !has_vcf {
            r.info("Suggested: Path D (ORA -> FASTQ -> BAM -> VCF)");
            r.say(&format!(
                "       First:  ./scripts/01-ora-to-fastq.sh {sm}   (run-all.sh does not read ORA)"
            ));
            r.say(&format!(
                "       Then:   ./scripts/run-all.sh {sm} <male|female>"
            ));
        } else if has_fastq && !has_bam && !has_vcf {
            r.info("Suggested: Path A (FASTQ -> BAM -> VCF)");
            r.say(&hint);
            r.say(&format!(
                "       or step by step, starting with:  ./scripts/02-alignment.sh {sm}"
            ));
        } else if has_bam && !has_vcf {
            r.info("Suggested: Path B (BAM -> VCF)");
            r.say(&hint);
            r.say(&format!(
                "       or step by step, starting with:  ./scripts/03-deepvariant.sh {sm}"
            ));
        } else if has_vcf {
            r.info("Suggested: Path C (VCF already available)");
            r.say(&hint);
            r.say(&format!(
                "       or step by step, starting with:  ./scripts/06-clinvar-screen.sh {sm}"
            ));
            if has_bam {
                r.say("       BAM also available — all pipeline steps can run.");
            } else if is_file(&bam) {
                r.say("       INPUT=vcf: the BAM is not read — BAM-dependent steps (4, 10, 15, 16, 18, 19, 20) will be skipped.");
            } else {
                r.say("       No BAM found — BAM-dependent steps (4, 10, 15, 16, 18, 19, 20) will be skipped.");
            }
        } else if !has_ora && !has_fastq && !has_bam && !has_vcf {
            if is_dir(&sd) {
                r.fail(&format!("No usable data found for sample '{sm}' in {sd}/"));
                r.say("       Expected one of:");
                r.say(&format!(
                    "         {sd}/fastq/{sm}_R1.fastq.gz + {sm}_R2.fastq.gz  (Path A)"
                ));
                r.say(&format!(
                    "         {sd}/aligned/{sm}_sorted.bam                         (Path B)"
                ));
                r.say(&format!(
                    "         {sd}/vcf/{sm}.vcf.gz                                 (Path C)"
                ));
                r.say(&format!(
                    "         {sd}/fastq/*.ora                                          (Path D)"
                ));
            } else {
                r.fail(&format!("Sample directory does not exist: {sd}/"));
                r.say("       Create it and place your data files inside:");
                r.say(&format!("         mkdir -p {sd}/fastq"));
                r.say("         # Copy your FASTQ/BAM/VCF files into the appropriate subdirectory");
            }
        }
    }

    /// `$(run_in SAMTOOLS samtools view -H BAM 2>/dev/null | grep "^@SQ" |
    /// grep PATTERN | head -1 | sed 's/.*LN://' | cut -f1 || echo "0")`
    /// under pipefail: the length after the last LN: of the first matching
    /// @SQ line, with "0" after it when a stage of the pipe failed.
    fn bam_chr1_len(&self, img: &str, bam: &str, pick: impl Fn(&str) -> bool) -> String {
        let (rc, out) = self.run_in(&[img, "samtools", "view", "-H", bam]);
        let first = lines(&out)
            .into_iter()
            .find(|l| l.starts_with("@SQ") && pick(l));
        let v = first.map(|l| {
            let after = match l.rfind("LN:") {
                Some(i) => &l[i + 3..],
                None => l,
            };
            after.split('\t').next().unwrap_or("").to_owned()
        });
        match (rc == 0, v) {
            (true, Some(v)) => v,
            (false, Some(v)) => format!("{v}\n0"),
            (_, None) => "0".to_owned(),
        }
    }
}

/// `grep -oE '[<,]length=[0-9]+' | cut -d= -f2` of one ##contig line: each
/// length, one per line.
fn contig_lengths(l: &str) -> String {
    let b = l.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if (b[i] == b'<' || b[i] == b',') && l[i + 1..].starts_with("length=") {
            let start = i + 1 + "length=".len();
            let end = start + l[start..].bytes().take_while(u8::is_ascii_digit).count();
            if end > start {
                out.push(&l[start..end]);
                i = end;
                continue;
            }
        }
        i += 1;
    }
    out.join("\n")
}

/// `find DIR -maxdepth 1 -name PATTERN | wc -l`: entries of DIR (any kind)
/// whose name matches.
fn count_in(dir: &str, m: impl Fn(&str) -> bool) -> usize {
    std::fs::read_dir(dir).map_or(0, |d| {
        d.flatten()
            .filter(|e| m(&e.file_name().to_string_lossy()))
            .count()
    })
}

/// The (pgs_id, licence) rows of assets/pgs_scores.tsv, as validate-setup.sh's
/// awk reads them: from the header line that starts with pgs_id on, every
/// line whose first column is PGS and digits.
pub fn pgs_rows(p: &Path) -> Vec<(String, String)> {
    let text = read_text(&p.to_string_lossy());
    let mut col: Option<usize> = None;
    let mut header = false;
    let mut out = Vec::new();
    for l in lines(&text) {
        let f: Vec<&str> = l.split('\t').collect();
        if !header && f[0] == "pgs_id" {
            col = f.iter().rposition(|c| *c == "licence");
            header = true;
            continue;
        }
        let id = f[0];
        if header
            && id.len() > 3
            && id.starts_with("PGS")
            && id[3..].bytes().all(|b| b.is_ascii_digit())
        {
            let lic = col.and_then(|c| f.get(c)).copied().unwrap_or("");
            out.push((id.to_owned(), lic.trim_end_matches('\t').to_owned()));
        }
    }
    out
}

/// `gzip -t FILE 2>/dev/null`, with the host's gzip as the shell ran it.
fn gzip_t(p: &str) -> bool {
    Command::new("gzip")
        .args(["-t", p])
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// consumer_sites_digest of common.sh: the sha256 of the table's content
/// read through gzip.
pub fn consumer_sites_digest(p: &str) -> String {
    let Ok(f) = std::fs::File::open(p) else {
        return String::new();
    };
    let mut d = flate2::read::MultiGzDecoder::new(std::io::BufReader::new(f));
    let mut h = <sha2::Sha256 as sha2::Digest>::new();
    let mut buf = vec![0u8; 1 << 16];
    loop {
        match d.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => sha2::Digest::update(&mut h, &buf[..n]),
            Err(_) => break,
        }
    }
    hex(&sha2::Digest::finalize(h))
}

/// `[[ X =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]`.
fn is_date(x: &str) -> bool {
    let b = x.as_bytes();
    b.len() == 10
        && b[4] == b'-'
        && b[7] == b'-'
        && b.iter()
            .enumerate()
            .all(|(i, c)| i == 4 || i == 7 || c.is_ascii_digit())
}

/// Whole days from DATE (YYYY-MM-DD, midnight UTC) to now, as `date -u -d`
/// reads it; `None` for a date that is not on the calendar.
pub fn days_since(date: &str) -> Option<i64> {
    let y: i64 = date[..4].parse().ok()?;
    let m: i64 = date[5..7].parse().ok()?;
    let d: i64 = date[8..].parse().ok()?;
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    let dim = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    if !(1..=12).contains(&m) || d < 1 || d > dim[(m - 1) as usize] {
        return None;
    }
    // Days from 1970-01-01 (Howard Hinnant's days_from_civil).
    let yy = if m <= 2 { y - 1 } else { y };
    let era = yy.div_euclid(400);
    let yoe = yy - era * 400;
    let mp = (m + 9) % 12;
    let doy = (153 * mp + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    let days = era * 146_097 + doe - 719_468;
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() as i64);
    Some((now - days * 86_400) / 86_400)
}

const HELP: &str = "usage: og setup check [SAMPLE]

  Checks the setup before a run, as scripts/validate-setup.sh did and
  still does through og: the system, GENOME_DIR, the reference data the
  component table names, the images of versions.env and their pinned tags,
  og's own binaries, and with SAMPLE the sample's reads, BAM or VCF (in the
  samtools and bcftools images). Exit status: 0 with no [FAIL], 1 with one
  or more, 2 for a sample name or a reference outside GENOME_DIR the scripts
  refuse. Reads GENOME_DIR, REF_FASTA, INPUT, TOOLS, ALLOW_ALT_REFERENCE,
  READ_COUNT_CHECK, NXF_WORK, OG_ARCHES, THREADS and CONTAINER_ENGINE.
";

/// `og setup check [SAMPLE]`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if matches!(args.first().map(String::as_str), Some("--help" | "-h")) {
        let _ = write!(out, "{HELP}");
        return 0;
    }
    let root = fetch::root_of(None);
    let sample = args.first().cloned().unwrap_or_default();
    let check = match Check::new(&root, Vars::from_process(), &sample) {
        Ok(c) => c,
        Err((code, m)) => {
            let _ = write!(err, "{m}");
            return code;
        }
    };
    if !sample.is_empty()
        && let Err(m) = validate_sample(&sample)
    {
        let _ = writeln!(err, "{m}");
        return 2;
    }
    let mut r = Report::new(out, Colors::for_stdout());
    let code = check.run(&mut r);
    r.flush();
    code
}

/// The installed paths of ROWS ((component key, its installed templates))
/// that break the rule of this module, each with what is wrong: a path no
/// check asked for (SEEN) and [`NOT_CHECKED`] does not list, a path both
/// read and listed, and a listed path no row has. Empty when every path is
/// read by a check or listed with its reason.
pub fn uncovered<'a>(
    rows: impl IntoIterator<Item = (&'a str, &'a [&'a str])>,
    seen: &BTreeSet<(String, String)>,
) -> Vec<String> {
    let listed: BTreeSet<(String, String)> = NOT_CHECKED
        .iter()
        .map(|(k, t, _)| ((*k).to_owned(), (*t).to_owned()))
        .collect();
    let mut all = BTreeSet::new();
    let mut out = Vec::new();
    for (key, templates) in rows {
        for t in templates {
            let p = (key.to_owned(), (*t).to_owned());
            match (seen.contains(&p), listed.contains(&p)) {
                (false, false) => out.push(format!(
                    "{key}: {t} is read by no check and not in NOT_CHECKED"
                )),
                (true, true) => out.push(format!(
                    "{key}: {t} is read by a check and also in NOT_CHECKED"
                )),
                _ => {}
            }
            all.insert(p);
        }
    }
    for (k, t) in listed.difference(&all) {
        out.push(format!("{k}: {t} is in NOT_CHECKED but in no row"));
    }
    for (k, t) in seen.difference(&all) {
        out.push(format!("{k}: {t} is read by a check but in no row"));
    }
    out
}
