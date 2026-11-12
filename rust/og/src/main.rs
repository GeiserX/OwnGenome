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

/// The release builds are static musl binaries; see `og/Cargo.toml`.
#[cfg(target_env = "musl")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn tree_hash() -> &'static str {
    &TREE_MARK["og-tree-hash:".len()..TREE_MARK.len() - 1]
}

const USAGE: &str = "usage: og <tool> [args]
       og map <verb> [args]
       og run state|tasks|eta|headroom|delta [args]
       og setup mask-fasta SRC BED DEST
       og setup unzip ZIP DEST
       og setup genetic-map --zip ZIP --out DIR [--rfmix] CHROM...
       og vcf csq-rank [options]
       og table constraint-join --gene-col NAME [--constrained] CONSTRAINT TABLE
       og sv svlen-fix|inv-filter [args]
       og pgx cpic-report|step36-outside [args]
       og pgx refblock-fill|mito-merge|cyp2d6-bed|cyp2d6-measure|cyp2d6-check [args]
       og clinvar join-adjacent|hits|deep-intronic [args]
       og str motif|catalog [args]
       og cdna check [args]
       og cdna scan [args]
       og comphet phase [args]
       og bam contigs BAM
       og qc sex-loss --table TABLE [--regions REGIONS] [--vcf VCF]
       og prs format [args]
       og chip prepare|report|check-table [args]
       og cn bed|call|intervals|check-bam|marker-regions|marker-read [args]
       og family plan|prepare|check|report|segments [args]
       og consumer build-sites|export|check [args]
       og ancestry-fine panel-map|check-map|order|sample-bed|chroms|report [args]
       og vntyper report [args]
       og locityper report [args]
       og mocha report [args]
       og mei report [args]
       og mei annotsv-input [args]
       og viral report [args]
       og y str-bed [args]
       og y str-report [args]
       og exomiser prepare [args]
       og exomiser report [args]
       og dev parity-compare|contract|schema|report-only|remnants [args]
       og version
       og --tree-hash
";

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
            Some("unzip") => ExitCode::from(og_setup::unzip::cli(&rest[1..])),
            Some("genetic-map") => ExitCode::from(og_setup::genetic_map::cli(&rest[1..])),
            Some(other) => usage_error(&format!("unknown setup tool '{other}'")),
            None => usage_error("the following arguments are required: setup tool"),
        },
        "qc" => match rest.first().map(String::as_str) {
            Some("sex-loss") => ExitCode::from(og_qc::sex_loss::cli(&rest[1..])),
            Some(other) => usage_error(&format!("unknown qc tool '{other}'")),
            None => usage_error("the following arguments are required: qc tool"),
        },
        "run" => ExitCode::from(og_runlog::cli::main(rest)),
        "consumer" => ExitCode::from(og_steps::consumer::cli(rest)),
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
                // cpic-report and step36-outside, and the usage for any other verb.
                None => ExitCode::from(og_pgx::main(rest)),
            },
            None => ExitCode::from(og_pgx::main(rest)),
        },
        "dev" => ExitCode::from(og_dev::main(rest)),
        "family" => ExitCode::from(og_steps::family::cli(rest)),
        "clinvar" => ExitCode::from(og_clinvar::cli::main(rest)),
        "chip" => match rest.first() {
            Some(verb) => ExitCode::from(og_chip::cli(
                verb,
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            None => usage_error("the following arguments are required: chip tool"),
        },
        "cn" => ExitCode::from(og_cn::cli::main(
            rest,
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
        "vntyper" | "locityper" | "mocha" | "mei" | "viral" | "y" | "exomiser" => {
            let verb = rest.first().map(String::as_str);
            let args = rest.get(1..).unwrap_or_default();
            let (out, err) = (&mut std::io::stdout().lock(), &mut std::io::stderr().lock());
            ExitCode::from(match tool.as_str() {
                "vntyper" => og_reports::vntyper::cli(verb, args, out, err),
                "locityper" => og_reports::locityper::cli(verb, args, out, err),
                "mocha" => og_reports::mocha::cli(verb, args, out, err),
                "mei" => og_reports::mei::cli(verb, args, out, err),
                "viral" => og_reports::viral::cli(verb, args, out, err),
                "y" => og_reports::ystr::cli(verb, args, out, err),
                _ => og_reports::exomiser::cli(verb, args, out, err),
            })
        }
        "ancestry-fine" => ExitCode::from(og_ancestry::cli(
            rest,
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
        "str" => match rest.first().map(String::as_str) {
            #[cfg(feature = "bam")]
            Some("motif") => ExitCode::from(og_bam::motif::cli(
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            #[cfg(not(feature = "bam"))]
            Some("motif") => no_bam("str motif"),
            Some("catalog") => ExitCode::from(og_setup::strchive::cli(&rest[1..])),
            _ => usage_error("og str takes motif or catalog"),
        },
        "prs" => match rest.first().map(String::as_str) {
            Some("format") => ExitCode::from(og_steps::prs::cli(&rest[1..])),
            _ => usage_error("og prs takes format"),
        },
        "cdna" | "comphet" | "bam" => {
            let verb = rest.first().map(String::as_str);
            let words = format!("{tool} {}", verb.unwrap_or_default());
            #[cfg(feature = "bam")]
            {
                let args = rest.get(1..).unwrap_or_default();
                let (out, err) = (&mut std::io::stdout().lock(), &mut std::io::stderr().lock());
                match words.as_str() {
                    "cdna check" => ExitCode::from(og_bam::cdna::check_cli(args, out, err)),
                    "cdna scan" => ExitCode::from(og_bam::cdna::scan_cli(args, out, err)),
                    "comphet phase" => ExitCode::from(og_bam::comphet::cli(args, out, err)),
                    "bam contigs" => ExitCode::from(og_bam::contigs_cli(args, out, err)),
                    _ => usage_error(&format!("unknown tool '{}'", words.trim_end())),
                }
            }
            #[cfg(not(feature = "bam"))]
            match words.as_str() {
                "cdna check" | "cdna scan" | "comphet phase" | "bam contigs" => no_bam(&words),
                _ => usage_error(&format!("unknown tool '{}'", words.trim_end())),
            }
        }
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
