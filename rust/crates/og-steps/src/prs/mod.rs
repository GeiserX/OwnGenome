//! `og prs format`: step 25's scoring files as pgsc_calc reads them
//! (`PRS_PREPARE`, `scripts/25-prs.sh`).
//!
//! ```text
//! og prs format --scores DIR [--ids A,B] [--labels TSV] --out DIR --alleles F
//!               [--keep-x] [--keep-ambiguous --reference FASTA]
//! ```
//!
//! Each GRCh38-harmonised PGS Catalog file of `--scores` (`*.txt.gz`,
//! `*.txt`) is written as `OUT/<pgs_id>.txt.gz`, pgsc_calc's custom GRCh38
//! file: `chr_name` and `chr_position` from `hm_chr` and `hm_pos`, the effect
//! and other allele, the weight, and the trait from the score list. Rows on
//! chromosomes 1 to 22 are kept, and chrX with `--keep-x`; chrY and MT never.
//! `OUT` also gets `chrx_rows.tsv`, `score_notes.tsv` and
//! `score_evaluations.tsv`, and `--alleles` every effect and other allele,
//! sorted and unique, for genotyping the score positions from the gVCF.
//!
//! With `--keep-ambiguous` a score keeps its strand-ambiguous rows (A/T, C/G)
//! when its other single-base rows agree with `--reference` on the forward
//! strand (at least one agrees, none only after a flip), and loses them when
//! they do not.
//!
//! Each score is streamed, so the memory is that of the allele set, not of
//! the scores. A refusal is `ERROR: ...` on stderr and exit 1; a bad command
//! line exits 2.

mod alleles;
mod labels;
mod score;

pub use labels::{Eval, Labels};
pub use score::{Ambiguous, Summary, pgs_id_of};

use alleles::Alleles;
use og_intervals::fasta::{FastaError, Indexed};
use og_io::AtomicFile;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

/// Beside the scores: per score its chrX rows and whether they were scored.
pub const CHRX_FILE: &str = "chrx_rows.tsv";
/// Beside the scores: per score the sex it applies to and its ambiguous rows.
pub const SCORE_NOTES_FILE: &str = "score_notes.tsv";
/// Beside the scores: per score its PGS Catalog evaluations.
pub const EVALS_FILE: &str = "score_evaluations.tsv";

const USAGE: &str = "usage: og prs format --scores DIR [--ids A,B] [--labels TSV] --out DIR --alleles F\n                     [--keep-x] [--keep-ambiguous --reference FASTA]\n";

/// The command line.
#[derive(Debug, Default)]
pub struct Args {
    pub scores: PathBuf,
    pub ids: String,
    pub labels: Option<PathBuf>,
    pub out: PathBuf,
    pub alleles: PathBuf,
    pub keep_x: bool,
    pub keep_ambiguous: bool,
    pub reference: Option<PathBuf>,
}

impl Args {
    /// argparse's reading: `--opt value` or `--opt=value`, the last one wins.
    pub fn parse(argv: &[String]) -> Result<Args, String> {
        let mut a = Args::default();
        let (mut scores, mut out, mut alleles) = (None, None, None);
        let mut it = argv.iter();
        while let Some(arg) = it.next() {
            let (key, inline) = match arg.split_once('=') {
                Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
                _ => (arg.as_str(), None),
            };
            let mut value = || {
                inline
                    .clone()
                    .or_else(|| it.next().cloned())
                    .ok_or_else(|| format!("argument {key}: expected one argument"))
            };
            match key {
                "--scores" => scores = Some(PathBuf::from(value()?)),
                "--ids" => a.ids = value()?,
                "--labels" => {
                    a.labels = Some(PathBuf::from(value()?)).filter(|p| !p.as_os_str().is_empty())
                }
                "--out" => out = Some(PathBuf::from(value()?)),
                "--alleles" => alleles = Some(PathBuf::from(value()?)),
                "--reference" => {
                    a.reference =
                        Some(PathBuf::from(value()?)).filter(|p| !p.as_os_str().is_empty())
                }
                "--keep-x" if inline.is_none() => a.keep_x = true,
                "--keep-ambiguous" if inline.is_none() => a.keep_ambiguous = true,
                _ => return Err(format!("unrecognized arguments: {arg}")),
            }
        }
        let missing: Vec<&str> = [
            ("--scores", scores.is_none()),
            ("--out", out.is_none()),
            ("--alleles", alleles.is_none()),
        ]
        .into_iter()
        .filter_map(|(k, m)| m.then_some(k))
        .collect();
        if !missing.is_empty() {
            return Err(format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ));
        }
        a.scores = scores.unwrap_or_default();
        a.out = out.unwrap_or_default();
        a.alleles = alleles.unwrap_or_default();
        Ok(a)
    }
}

/// `og prs format ARGS`: the exit code.
pub fn cli(argv: &[String]) -> u8 {
    let args = match Args::parse(argv) {
        Ok(a) => a,
        Err(e) => {
            eprint!("{USAGE}");
            eprintln!("og prs format: error: {e}");
            return 2;
        }
    };
    let stdout = io::stdout();
    match run(&args, &mut stdout.lock()) {
        Ok(()) => 0,
        Err(e) => {
            let _ = stdout.lock().flush();
            eprintln!("ERROR: {e}");
            1
        }
    }
}

/// The scoring files of `dir`, sorted: `*.txt.gz` and `*.txt`, hidden ones
/// left out, as `sorted(glob(...))` lists them. A folder that cannot be
/// read has none.
fn score_files(dir: &Path) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter(|e| {
            let n = e.file_name();
            let n = n.to_string_lossy();
            !n.starts_with('.') && (n.ends_with(".txt.gz") || n.ends_with(".txt"))
        })
        .map(|e| dir.join(e.file_name()))
        .collect();
    files.sort();
    files
}

/// The file name, for messages.
fn base(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

fn two_files(pid: &str, dir: &Path, a: &Path, b: &Path) -> String {
    let (a, b) = (base(a), base(b));
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    format!("{pid} has two scoring files in {}: {a}, {b}", dir.display())
}

/// The run; the error is the message after `ERROR: `.
pub fn run(a: &Args, out: &mut impl Write) -> Result<(), String> {
    if a.keep_ambiguous && a.reference.is_none() {
        return Err("--keep-ambiguous needs --reference, the FASTA the strand check reads".into());
    }
    let labels = match &a.labels {
        Some(p) => Labels::read(p)?,
        None => Labels::default(),
    };
    let mut fasta = match (&a.reference, a.keep_ambiguous) {
        (Some(r), true) => Some(Indexed::open(r).map_err(|e| match e {
            FastaError::Compressed(_) => format!(
                "{} is compressed; the strand check reads an uncompressed FASTA",
                r.display()
            ),
            e => format!(
                "cannot read the reference {} and its .fai: {e}",
                r.display()
            ),
        })?),
        _ => None,
    };

    let all = score_files(&a.scores);
    let mut by_id: HashMap<String, Vec<&PathBuf>> = HashMap::new();
    for p in &all {
        by_id.entry(pgs_id_of(p)).or_default().push(p);
    }
    let files: Vec<&PathBuf> = if a.ids.is_empty() {
        all.iter().collect()
    } else {
        let want: Vec<&str> = a
            .ids
            .split(',')
            .map(|x| pycompat::text::strip(x, None))
            .filter(|x| !x.is_empty())
            .collect();
        let missing: Vec<&str> = want
            .iter()
            .copied()
            .filter(|x| !by_id.contains_key(*x))
            .collect();
        if !missing.is_empty() {
            return Err(format!(
                "no scoring file for {} in {}",
                missing.join(", "),
                a.scores.display()
            ));
        }
        want.iter().map(|x| by_id[*x][0]).collect()
    };
    for p in &files {
        if let [x, y, ..] = by_id[&pgs_id_of(p)].as_slice() {
            return Err(two_files(&pgs_id_of(p), &a.scores, x, y));
        }
    }
    if files.is_empty() {
        return Err(format!(
            "no scoring file (*.txt.gz or *.txt) in {}",
            a.scores.display()
        ));
    }
    fs::create_dir_all(&a.out).map_err(|e| format!("{}: {e}", a.out.display()))?;

    let label = |pid: &str| labels.trait_of(pid).map(str::to_owned);
    let opts = score::Options {
        keep_x: a.keep_x,
        label: &label,
    };
    let mut set = Alleles::default();
    let mut done: Vec<(Summary, &PathBuf)> = Vec::new();
    for p in files {
        let s = score::format(p, &a.out, &opts, fasta.as_mut(), &mut set)?;
        if let Some((_, q)) = done.iter().find(|(d, q)| d.pid == s.pid && *q != p) {
            return Err(two_files(&s.pid, &a.scores, q, p));
        }
        let kept = if a.keep_x && s.chrx > 0 {
            format!(", {} on chrX scored", s.chrx)
        } else {
            String::new()
        };
        let dropped = if s.off_autosomes > 0 {
            format!(", {} off the autosomes dropped", s.off_autosomes)
        } else {
            String::new()
        };
        let mut ambig = if s.ambiguous > 0 {
            format!(", {} strand-ambiguous {}", s.ambiguous, s.fate.as_str())
        } else {
            String::new()
        };
        if fasta.is_some() {
            ambig += &format!(
                " (strand check: {} rows agree with the reference, {} flipped)",
                s.agree, s.flipped
            );
        }
        let sex = labels
            .applies_to(&s.pid)
            .map_or_else(String::new, |x| format!(", {x} only"));
        let label = labels.trait_of(&s.pid).unwrap_or("label from the file");
        writeln!(
            out,
            "  {}: {} GRCh38 rows{kept}{dropped}{ambig}{sex} ({label})",
            s.pid, s.written
        )
        .map_err(|e| e.to_string())?;
        done.push((s, p));
    }

    let write = |name: &str, body: String| -> Result<(), String> {
        let p = a.out.join(name);
        fs::write(&p, body).map_err(|e| format!("{}: {e}", p.display()))
    };
    let chrx = if a.keep_x { "scored" } else { "left out" };
    let mut text = String::from("pgs_id\tchrx_rows\tchrx\n");
    for (s, _) in &done {
        text += &format!("{}\t{}\t{chrx}\n", s.pid, s.chrx);
    }
    write(CHRX_FILE, text)?;
    let mut text = String::from("pgs_id\tapplies_to\tambiguous_rows\tambiguous\n");
    for (s, _) in &done {
        text += &format!(
            "{}\t{}\t{}\t{}\n",
            s.pid,
            labels.applies_to(&s.pid).unwrap_or(""),
            s.ambiguous,
            s.fate.as_str()
        );
    }
    write(SCORE_NOTES_FILE, text)?;
    let mut text = String::from("pgs_id\tppm_id\tancestry\teffect\n");
    for (s, _) in &done {
        let ev = labels.evals(&s.pid);
        if !ev.is_empty() {
            let join = |f: fn(&Eval) -> &str| ev.iter().map(f).collect::<Vec<_>>().join(";");
            text += &format!(
                "{}\t{}\t{}\t{}\n",
                s.pid,
                join(|e| &e.ppm),
                join(|e| &e.ancestry),
                join(|e| &e.effect)
            );
        }
    }
    write(EVALS_FILE, text)?;

    let fail = |e: io::Error| format!("{}: {e}", a.alleles.display());
    let mut f = AtomicFile::create(&a.alleles).map_err(fail)?;
    set.write(&mut f).map_err(fail)?;
    f.commit().map_err(fail)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// A fresh folder under the system temp folder, removed on drop.
    struct Dir(PathBuf);

    impl Dir {
        fn new(tag: &str) -> Dir {
            let d = std::env::temp_dir().join(format!("og-prs-{tag}-{}", std::process::id()));
            let _ = fs::remove_dir_all(&d);
            fs::create_dir_all(&d).expect("dir");
            Dir(d)
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn gz(path: &Path, text: &str) {
        let mut w = flate2::write::GzEncoder::new(
            fs::File::create(path).expect("create"),
            flate2::Compression::fast(),
        );
        w.write_all(text.as_bytes()).expect("write");
        w.finish().expect("finish");
    }

    fn ungz(path: &Path) -> String {
        let mut s = String::new();
        flate2::read::MultiGzDecoder::new(fs::File::open(path).expect("open"))
            .read_to_string(&mut s)
            .expect("read");
        s
    }

    fn args(d: &Path, extra: &[&str]) -> Args {
        let mut v: Vec<String> = ["--scores", "in", "--out", "pgs", "--alleles", "alleles.tsv"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        v.extend(extra.iter().map(|s| (*s).to_owned()));
        let mut a = Args::parse(&v).expect("args");
        a.scores = d.join(&a.scores);
        a.out = d.join(&a.out);
        a.alleles = d.join(&a.alleles);
        a.labels = a.labels.map(|p| d.join(p));
        a.reference = a.reference.map(|p| d.join(p));
        a
    }

    #[test]
    fn writes_the_score_its_notes_and_the_sorted_alleles() {
        let d = Dir::new("write");
        fs::create_dir_all(d.0.join("in")).expect("in");
        // The header goes on after the column line, as the Python allowed.
        gz(
            &d.0.join("in/PGS1_hmPOS_GRCh38.txt.gz"),
            "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\thm_inferOtherAllele\n\
             chr10\t5\ta\tA/G\t0.5\tg\n2\t7\tC\tT\t-1\t\nX\t9\tA\tG\t1\t\nY\t1\tA\tG\t1\t\n1\t3\tA\tT\t2\t\n\
             #pgs_id=PGS1\n#HmPOS_build=GRCh38\n",
        );
        fs::write(
            d.0.join("labels.tsv"),
            "pgs_id\ttrait_reported\tapplies_to\nPGS1\tA = trait\tFemale\n",
        )
        .expect("labels");
        let mut out = Vec::new();
        run(
            &args(&d.0, &["--labels", "labels.tsv", "--keep-x"]),
            &mut out,
        )
        .expect("run");
        assert_eq!(
            String::from_utf8(out).expect("utf8"),
            "  PGS1: 4 GRCh38 rows, 1 on chrX scored, 1 off the autosomes dropped, 1 strand-ambiguous dropped, female only (A = trait)\n"
        );
        assert_eq!(
            ungz(&d.0.join("pgs/PGS1.txt.gz")),
            "#pgs_id=PGS1\n#pgs_name=PGS1\n#trait_reported=A - trait\n#genome_build=GRCh38\n\
             chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n\
             10\t5\tA\tG\t0.5\n2\t7\tC\tT\t-1\nX\t9\tA\tG\t1\n1\t3\tA\tT\t2\n"
        );
        assert_eq!(
            fs::read_to_string(d.0.join("alleles.tsv")).expect("alleles"),
            "chr1\t3\tA\nchr1\t3\tT\nchr10\t5\tA\nchr10\t5\tG\nchr2\t7\tC\nchr2\t7\tT\nchrX\t9\tA\nchrX\t9\tG\n"
        );
        assert_eq!(
            fs::read_to_string(d.0.join("pgs/score_notes.tsv")).expect("notes"),
            "pgs_id\tapplies_to\tambiguous_rows\tambiguous\nPGS1\tfemale\t1\tdropped\n"
        );
        assert_eq!(
            fs::read_to_string(d.0.join("pgs/chrx_rows.tsv")).expect("chrx"),
            "pgs_id\tchrx_rows\tchrx\nPGS1\t1\tscored\n"
        );
    }

    #[test]
    fn refuses_two_files_of_one_id_before_writing() {
        let d = Dir::new("two");
        fs::create_dir_all(d.0.join("in")).expect("in");
        let score =
            "#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n1\t5\tA\t1\n";
        fs::write(d.0.join("in/PGS2.txt"), score).expect("txt");
        gz(&d.0.join("in/PGS2.txt.gz"), score);
        let err = run(&args(&d.0, &[]), &mut Vec::new()).expect_err("refused");
        assert!(
            err.starts_with("PGS2 has two scoring files in ")
                && err.ends_with(": PGS2.txt, PGS2.txt.gz"),
            "{err}"
        );
        assert!(!d.0.join("pgs").exists());
    }

    #[test]
    fn strand_check_against_a_reference() {
        let d = Dir::new("strand");
        fs::create_dir_all(d.0.join("in")).expect("in");
        // chr1 is GAAGccAT: soft-masked bases count, upper-cased.
        fs::write(d.0.join("ref.fa"), ">chr1\nGAAG\nccAT\n").expect("fa");
        fs::write(d.0.join("ref.fa.fai"), "chr1\t8\t6\t4\t5\n").expect("fai");
        let head = "#pgs_id=PGS3\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n";
        let rows = "1\t2\tA\tG\t1\n1\t3\tA\tT\t1\n1\t6\tT\tC\t1\n";
        gz(&d.0.join("in/PGS3.txt.gz"), &format!("{head}{rows}"));
        let mut out = Vec::new();
        run(
            &args(&d.0, &["--keep-ambiguous", "--reference", "ref.fa"]),
            &mut out,
        )
        .expect("run");
        assert!(String::from_utf8(out).expect("utf8").contains(
            "1 strand-ambiguous kept (strand check: 2 rows agree with the reference, 0 flipped)"
        ));
        // chr1:7 is A: T/C agrees only on the other strand.
        gz(
            &d.0.join("in/PGS3.txt.gz"),
            &format!("{head}{rows}1\t7\tT\tC\t1\n"),
        );
        run(
            &args(&d.0, &["--keep-ambiguous", "--reference", "ref.fa"]),
            &mut Vec::new(),
        )
        .expect("run");
        assert!(!ungz(&d.0.join("pgs/PGS3.txt.gz")).contains("\t3\tA\tT\t"));
        assert!(
            fs::read_to_string(d.0.join("pgs/score_notes.tsv"))
                .expect("notes")
                .ends_with("\t1\tdropped, strand check failed\n")
        );
    }

    #[test]
    fn command_line_errors() {
        let e =
            |v: &[&str]| Args::parse(&v.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>()).err();
        assert_eq!(
            e(&["--scores", "x"]).as_deref(),
            Some("the following arguments are required: --out, --alleles")
        );
        assert_eq!(
            e(&["--scores"]).as_deref(),
            Some("argument --scores: expected one argument")
        );
        assert_eq!(
            e(&["--keep-x=1"]).as_deref(),
            Some("unrecognized arguments: --keep-x=1")
        );
        assert!(e(&["--scores=a", "--out", "b", "--alleles", "c"]).is_none());
    }
}
