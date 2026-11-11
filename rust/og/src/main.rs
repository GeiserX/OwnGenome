//! `og`: OwnGenome's one binary. `og <tool> [args]`; each ported Python
//! tool becomes a subcommand named after its file in `bin/`.

use std::io::Write;
use std::process::ExitCode;

mod map;
mod probe;

/// The tree hash between a fixed marker and a `;`, so a host check can read
/// it from the file without running it (`grep -ao 'og-tree-hash:[0-9a-f]*;'`):
/// an Apple Silicon host cannot run the x86_64 Linux binary. The `;` ends
/// it where the next string in the binary would otherwise run on.
const TREE_MARK: &str = concat!("og-tree-hash:", env!("OG_TREE_HASH"), ";");
const TARGET: &str = env!("OG_TARGET");

fn tree_hash() -> &'static str {
    &TREE_MARK["og-tree-hash:".len()..TREE_MARK.len() - 1]
}

const USAGE: &str = "usage: og <tool> [args]\n       og map <verb> [args]\n       og run state|tasks|eta|headroom|delta [args]\n       og setup mask-fasta SRC BED DEST\n       og vcf csq-rank [options]\n       og table constraint-join --gene-col NAME [--constrained] CONSTRAINT TABLE\n       og sv svlen-fix|inv-filter [args]\n       og pgx refblock-fill|mito-merge|cyp2d6-bed|cyp2d6-measure|cyp2d6-check [args]\n       og str motif [args]\n       og chip prepare|report|check-table [args]\n       og family plan|prepare|check|report|segments [args]\n       og dev parity-compare KIND ORACLE NEW\n       og dev remnants [--list FILE] [--root DIR]\n       og version\n       og --tree-hash\n";

/// argparse's way to refuse a command line: usage and the error on stderr, exit 2.
fn usage_error(msg: &str) -> ExitCode {
    eprint!("{USAGE}");
    eprintln!("og: error: {msg}");
    ExitCode::from(2)
}

/// A BAM tool called on a build without the `bam` feature (the macOS build
/// for a live host, rust/og.sh): it says so and exits 2, before any input
/// is read, so a caller never mistakes it for a result.
#[cfg(not(feature = "bam"))]
fn no_bam(tool: &str) -> ExitCode {
    eprintln!(
        "og: error: {tool} reads BAMs, and this og ({TARGET}) is built without BAM support; run it in the pipeline's Linux containers"
    );
    ExitCode::from(2)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(tool) = args.first() else {
        return usage_error("the following arguments are required: tool");
    };
    let rest = &args[1..];
    match tool.as_str() {
        "-h" | "--help" => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        "version" | "--version" if rest.is_empty() => {
            println!("og rust tree {} ({TARGET})", tree_hash());
            ExitCode::SUCCESS
        }
        "--tree-hash" if rest.is_empty() => {
            println!("{}", tree_hash());
            ExitCode::SUCCESS
        }
        // Hidden: the differential test against CPython drives it.
        "pycompat-probe" if rest.is_empty() => {
            let stdin = std::io::stdin().lock();
            let stdout = std::io::stdout().lock();
            match probe::run(stdin, stdout) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    let _ = writeln!(std::io::stderr(), "og pycompat-probe: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        "map" => map::run(rest, tree_hash()),
        "setup" => match rest.first().map(String::as_str) {
            Some("mask-fasta") => ExitCode::from(og_intervals::mask::cli(&rest[1..])),
            Some(other) => usage_error(&format!("unknown setup tool '{other}'")),
            None => usage_error("the following arguments are required: setup tool"),
        },
        "run" => ExitCode::from(og_runlog::cli::main(rest)),
        "vcf" => match rest.first().map(String::as_str) {
            Some("csq-rank") => ExitCode::from(og_vcf::rank::cli(&rest[1..])),
            Some(other) => usage_error(&format!("unknown vcf tool '{other}'")),
            None => usage_error("the following arguments are required: vcf tool"),
        },
        "table" => match rest.first().map(String::as_str) {
            Some("constraint-join") => ExitCode::from(og_steps::constraint::cli(&rest[1..])),
            Some(other) => usage_error(&format!("unknown table tool '{other}'")),
            None => usage_error("the following arguments are required: table tool"),
        },
        "sv" => match rest.first().map(String::as_str) {
            Some("svlen-fix") => ExitCode::from(og_filters::svlen::cli(&rest[1..])),
            Some("inv-filter") => ExitCode::from(og_filters::inv::cli(&rest[1..])),
            Some(other) => usage_error(&format!("unknown sv tool '{other}'")),
            None => usage_error("the following arguments are required: sv tool"),
        },
        "pgx" => match rest.first().map(String::as_str) {
            Some("refblock-fill") => ExitCode::from(og_filters::refblock::cli(&rest[1..])),
            Some("mito-merge") => ExitCode::from(og_filters::mito::cli(&rest[1..])),
            Some(other) => match other.strip_prefix("cyp2d6-") {
                #[cfg(feature = "bam")]
                Some(verb) => ExitCode::from(og_bam::cyp2d6::cli(
                    verb,
                    &rest[1..],
                    &mut std::io::stdout().lock(),
                    &mut std::io::stderr().lock(),
                )),
                #[cfg(not(feature = "bam"))]
                Some(_) => no_bam(&format!("pgx {other}")),
                None => usage_error(&format!("unknown pgx tool '{other}'")),
            },
            None => usage_error("the following arguments are required: pgx tool"),
        },
        "dev" => ExitCode::from(og_dev::main(rest)),
        "family" => ExitCode::from(og_steps::family::cli(rest)),
        "chip" => match rest.first() {
            Some(verb) => ExitCode::from(og_chip::cli(
                verb,
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            None => usage_error("the following arguments are required: chip tool"),
        },
        "str" => match rest.first().map(String::as_str) {
            #[cfg(feature = "bam")]
            Some("motif") => ExitCode::from(og_bam::motif::cli(
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            #[cfg(not(feature = "bam"))]
            Some("motif") => no_bam("str motif"),
            _ => usage_error("og str takes motif"),
        },
        // Hidden: the fetch smoke compares it with pysam (crates/og-bam).
        #[cfg(feature = "bam")]
        "bam-fetch-probe" if rest.len() == 2 => {
            let run = || -> std::io::Result<()> {
                let regions = std::io::BufReader::new(std::fs::File::open(&rest[1])?);
                og_bam::fetch_probe(rest[0].as_ref(), regions, std::io::stdout().lock())
            };
            match run() {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    let _ = writeln!(std::io::stderr(), "og bam-fetch-probe: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        #[cfg(not(feature = "bam"))]
        "bam-fetch-probe" => no_bam("bam-fetch-probe"),
        #[cfg(feature = "bam")]
        "bam-fetch-probe" => usage_error("bam-fetch-probe takes BAM REGIONS"),
        "version" | "--version" | "--tree-hash" | "pycompat-probe" => {
            usage_error(&format!("unrecognized arguments: {}", rest.join(" ")))
        }
        other => usage_error(&format!("unknown tool '{other}'")),
    }
}
