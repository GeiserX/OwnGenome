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

const USAGE: &str = "usage: og <tool> [verb] [args]
       og <tool> [verb] --help      one command's usage, on stdout

Run a genome (docs/getting-started.md):
  GENOME_DIR=/data ./scripts/run-all.sh SAMPLE male|female [nextflow options]
                       start a run; run it again to resume, it reuses every finished step
  ./scripts/validate-setup.sh [SAMPLE]
                       check the setup before a run
  og run state|tasks|eta|headroom|delta [args]
                       inspect a run's progress and time left; og run never starts a run
  og recheck --sample S [--genome-dir DIR] [--offline FILE] [--rerun [--refresh-clinvar]] | --affected KEYS | --launch-dir DIR
                       say which results of a finished run a newer ClinVar or catalog would change
  og lookup query|resolve|report [args]
                       your genotype at an rsID, a position or in a gene, from a finished run
  og report render [args]
                       the text and HTML report from a run's summary.json
  og live serve [args]
                       the live site of a run
  og view serve [args]
                       igv.js over one sample's reads and calls

Pipeline helpers (the steps run these):
  og map build --root DIR --run RUN --assets DIR --cache-dir DIR [--predicate-table]
  og map genes --root DIR --cache-dir DIR (--chrom CHROM | --names)
  og map svg [args]
  og setup mask-fasta SRC BED DEST
  og setup unzip ZIP DEST
  og setup genetic-map --zip ZIP --out DIR CHROM...
  og setup chrm-shift [args]
  og setup revel [args]
  og setup roh-af [args]
  og vcf csq-rank [options]
  og vcf index [args]
  og vcf shift-back [args]
  og table constraint-join --gene-col NAME [--constrained] CONSTRAINT TABLE
  og sv svlen-fix|inv-filter [args]
  og sv add-callers [args]
  og sv postfilter [args]
  og sv survivor-prep [args]
  og sv consensus-check [args]
  og pgx cpic-report|step36-outside [args]
  og pgx refblock-fill|mito-merge|cyp2d6-bed|cyp2d6-measure|cyp2d6-check [args]
  og pgx outside-calls|hla-concordance|pypgx-summary [args]
  og paraphase check-bam|table [args]
  og clinvar join-adjacent|hits|deep-intronic|conflicting-subset|conflicting-hits [args]
  og str motif|catalog [args]
  og cdna check [args]
  og cdna scan [args]
  og comphet phase [args]
  og bam contigs BAM
  og qc sex-loss --table TABLE [--regions REGIONS] [--vcf VCF]
  og qc sample [args]
  og prs format [args]
  og prs table [args]
  og mito annotate|deletions|copies [args]
  og summary collect [--only {a,b,c}] [args]
  og report index [args]
  og report charts [args]
  og live results [args]
  og live family [args]
  og live demo [args]
  og couple screen [args]
  og couple genes [args]
  og facts --sample S --sample-dir DIR --reference FASTA [--assets DIR] [--out FILE]
  og roh summary --sample SAMPLE --roh ROH [--out OUT]
  og roh af-file --panel PANEL --group GROUP --out OUT
  og chip prepare|report|check-table [args]
  og cn bed|call|intervals|check-bam|marker-regions|marker-read [args]
  og traits table [args]
  og traits phase-sites [args]
  og traits check [args]
  og bloodtype report|check-markers|systems [args]
  og family plan|prepare|check|report|segments [args]
  og family ibd [args]
  og consumer build-sites|export|check [args]
  og ancestry-fine panel-map|check-map|chroms|report [args]
  og ancestry-fine freqs [args]
  og ancestry-fine project [args]
  og ancestry-fine paint [args]
  og ancestry-fine simulate [args]
  og archaic table|strict|genotypes|report [args]
  og archaic compare [args]
  og archaic genomes [args]
  og ancient select|lift|person-geno|align|report [args]
  og vntyper report [args]
  og locityper report [args]
  og mocha report [args]
  og mei report [args]
  og mei annotsv-input [args]
  og viral report [args]
  og y str-bed [args]
  og y str-report [args]
  og y tree-positions --tree-dir DIR OUT
  og y tree-place [args]
  og y private-snps [args]
  og y yleaf-positions --data DATA OUT
  og exomiser prepare [args]
  og exomiser report [args]
  og longread clair3-model [args]

This og:
  og version
  og --tree-hash

Developer and CI tools: og dev --help lists them.
";

/// The developer and CI tools, which `og dev --help` prints. They stay out of
/// [`USAGE`], so `og --help` lists only what a user or a step runs.
const DEV_USAGE: &str = "  og dev parity-compare|contract|schema|report-only|remnants [args]
  og dev scan personal-data|images [args]
  og dev versions check-run|expand|coupled [args]
  og dev pgs trait|label-match|metric-match [args]
  og dev fixture comphet OUT [--swap]
  og dev fixture trim-pharmcat REPORT OUT [MAX]
  og dev evidence clashes RAW.vcf.gz SPLIT.vcf.gz DATA_DIR OUT.tsv
  og dev evidence hg002-gvcf TRUTH.vcf.gz TRUTH.bed REF.fasta SITES.tsv OUT.vcf
  og dev evidence y-bench TREE.json CALLS.tsv OUT.tsv MIN
  og dev renovate-check LOG
  og dev pharmcat diplotypes REPORT [--called|--called-lines]
  og dev nf-trace requests|attempts TRACE [options]
  og dev smoke mini-ref|long-reads|eh-catalog|csv-cell|check [args]
  og dev helpers [args]
  og dev build cn-intervals|callability|sv-rule-regions [args]
  og dev build roh-af [args]
  og dev build acmg-lof --clingen FILE --overrides FILE --ttn-exons FILE --table FILE [--check]
  og dev demo-sample [args]
  og dev demo-events [args]
  og dev params-doc [args]
  og dev python-gate [args]
  og dev everything-check [args]
  og dev licences [args]
  og dev docs-check [args]
  og dev stub-http [args]
  og dev steps-check [args]
  og dev traits-sources [args]
  og dev pictures take|live|check|same-sample [args]
  og dev freshness [--out FILE] [--sections a,b] | --update-issue BODY --repo O/N | --self-test
";

/// Commands whose own parser answers `-h`/`--help` with a fuller usage than
/// their line in [`USAGE`] or [`DEV_USAGE`]: `tool` for every verb of it,
/// `tool a|b` for those verbs. Any other command's help request is answered
/// from its lines there by [`help`], so a new verb answers `--help` with no
/// change here.
const OWN_HELP: &[&str] = &[
    "run",
    "recheck",
    "clinvar",
    "lookup",
    "family",
    "chip",
    "traits",
    "bloodtype",
    "cn",
    "paraphase",
    "ancient",
    "pgx cpic-report|step36-outside|refblock-fill|mito-merge|outside-calls|hla-concordance|pypgx-summary",
    "vcf csq-rank|index",
    "table constraint-join",
    "sv svlen-fix|inv-filter|add-callers|postfilter",
    "str catalog",
    "qc sex-loss",
    "report render",
    "couple screen|genes",
    "ancestry-fine panel-map|check-map|freqs|project|chroms|report|simulate",
    "archaic table|strict|genotypes|report",
    "vntyper report",
    "locityper report",
    "mocha report",
    "mei report|annotsv-input",
    "viral report",
    "y str-bed|str-report|tree-positions|tree-place|private-snps|yleaf-positions",
    "exomiser prepare|report",
    "dev freshness|contract|report-only|remnants|renovate-check|helpers|params-doc|traits-sources",
];

/// `og WORDS... -h|--help [...]`: the usage lines of the command named by
/// the words before the first option, on stdout, exit 0; the lines of its
/// verb when one matches, else every line of the tool. `None` leaves the
/// request to the command: it is in [`OWN_HELP`] (its `-h` then reaches it
/// as `--help`, which every one of them reads), names no tool with a line
/// (an unknown tool keeps its refusal), or holds no help word before `--`.
fn help(raw: &mut [String]) -> Option<ExitCode> {
    let end = raw.iter().position(|a| a == "--").unwrap_or(raw.len());
    let args: Vec<&str> = raw[..end].iter().map(String::as_str).collect();
    if !args.iter().any(|a| *a == "-h" || *a == "--help") {
        return None;
    }
    // A root option with its own usage line (`og --tree-hash`) is a tool here.
    let words: Vec<&str> = args
        .iter()
        .copied()
        .enumerate()
        .take_while(|(i, a)| !a.starts_with('-') || (*i == 0 && *a == "--tree-hash"))
        .map(|(_, a)| a)
        .collect();
    let (tool, verb) = (*words.first()?, words.get(1).copied());
    let own = OWN_HELP.iter().any(|e| match e.split_once(' ') {
        None => *e == tool,
        Some((t, vs)) => t == tool && verb.is_some_and(|v| vs.split('|').any(|x| x == v)),
    });
    if own {
        for a in &mut raw[..end] {
            if a == "-h" {
                *a = "--help".to_owned();
            }
        }
        return None;
    }
    let lines: Vec<(&str, Vec<&str>)> = USAGE
        .lines()
        .chain(DEV_USAGE.lines())
        .map(|l| l.trim_start_matches("usage:").trim())
        .filter_map(|l| {
            let mut w = l.split(' ');
            (w.next() == Some("og") && w.next() == Some(tool)).then(|| {
                (
                    l,
                    w.next().map_or_else(Vec::new, |v| v.split('|').collect()),
                )
            })
        })
        .collect();
    let of_verb: Vec<&str> = lines
        .iter()
        .filter(|(_, vs)| verb.is_some_and(|v| vs.contains(&v)))
        .map(|(l, _)| *l)
        .collect();
    let show = if of_verb.is_empty() {
        lines.iter().map(|(l, _)| *l).collect()
    } else {
        of_verb
    };
    if show.is_empty() {
        return None;
    }
    for (i, l) in show.iter().enumerate() {
        println!("{}{l}", if i == 0 { "usage: " } else { "       " });
    }
    Some(ExitCode::SUCCESS)
}

/// argparse's way to refuse a command line: usage and the error on stderr, exit 2.
fn usage_error(msg: &str) -> ExitCode {
    eprint!("{USAGE}");
    eprintln!("og: error: {msg}");
    ExitCode::from(2)
}

/// Runs a verb that writes to stdout and stderr; returns its exit code.
fn io(f: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8, args: &[String]) -> ExitCode {
    ExitCode::from(f(
        args,
        &mut std::io::stdout().lock(),
        &mut std::io::stderr().lock(),
    ))
}

/// Runs a server verb. Its request threads print to stderr themselves (the
/// access log, a summary.json it cannot read), so no lock on stdout or
/// stderr is held for the whole command: a held lock would block each of
/// those prints, and the request with it, for good.
fn serve(f: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8, args: &[String]) -> ExitCode {
    ExitCode::from(f(args, &mut std::io::stdout(), &mut std::io::stderr()))
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
    // og is built with panic = "unwind" for the sections of og report
    // render; a panic anywhere else still ends og at once, as abort did, so
    // no worker thread can die quietly and leave a cut output.
    og_render::chapters::abort_outside_guarded();
    // Refuse raw-byte scanner arguments before std::env::args can panic or
    // any path can be replaced with a different name. Other verbs retain
    // their existing argument handling.
    let raw: Vec<_> = std::env::args_os().skip(1).collect();
    if raw
        .get(..3)
        .is_some_and(|prefix| prefix == ["dev", "scan", "personal-data"])
        && raw.iter().any(|s| s.to_str().is_none())
    {
        eprintln!("ERROR: path is not valid UTF-8");
        return ExitCode::from(2);
    }
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    if let Some(code) = help(&mut args) {
        return code;
    }
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
        "map" if rest.first().is_some_and(|v| v == "svg") => io(og_map::svg::cli, &rest[1..]),
        "map" => map::run(rest, tree_hash()),
        "setup" => match rest.first().map(String::as_str) {
            Some("mask-fasta") => ExitCode::from(og_intervals::mask::cli(&rest[1..])),
            Some("unzip") => ExitCode::from(og_setup::unzip::cli(&rest[1..])),
            Some("genetic-map") => ExitCode::from(og_setup::genetic_map::cli(&rest[1..])),
            Some("chrm-shift") => io(og_vcf::shift::chrm_shift_cli, &rest[1..]),
            Some("revel") => io(og_setup::revel::cli, &rest[1..]),
            Some("roh-af") => io(og_setup::roh_af::cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown setup tool '{other}'")),
            None => usage_error("the following arguments are required: setup tool"),
        },
        "qc" => match rest.first().map(String::as_str) {
            Some("sex-loss") => ExitCode::from(og_qc::sex_loss::cli(&rest[1..])),
            Some("sample") => ExitCode::from(og_summary::qc::cli(
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            Some(other) => usage_error(&format!("unknown qc tool '{other}'")),
            None => usage_error("the following arguments are required: qc tool"),
        },
        "run" => ExitCode::from(og_runlog::cli::main(rest)),
        "consumer" => ExitCode::from(og_steps::consumer::cli(rest)),
        "vcf" => match rest.first().map(String::as_str) {
            Some("csq-rank") => ExitCode::from(og_vcf::rank::cli(&rest[1..])),
            Some("index") => io(og_io::tabix_write::cli, &rest[1..]),
            Some("shift-back") => io(og_vcf::shift::shift_back_cli, &rest[1..]),
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
            Some("add-callers") => ExitCode::from(og_steps::sv::add_callers::cli(&rest[1..])),
            Some("postfilter") => ExitCode::from(og_steps::sv::postfilter::cli(&rest[1..])),
            Some("survivor-prep") => io(og_steps::sv::survivor_prep::cli, &rest[1..]),
            Some("consensus-check") => io(og_steps::sv::consensus_check::cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown sv tool '{other}'")),
            None => usage_error("the following arguments are required: sv tool"),
        },
        "pgx" => match rest.first().map(String::as_str) {
            Some("refblock-fill") => ExitCode::from(og_filters::refblock::cli(&rest[1..])),
            Some("mito-merge") => ExitCode::from(og_filters::mito::cli(&rest[1..])),
            Some(verb @ ("outside-calls" | "hla-concordance" | "pypgx-summary")) => {
                let (out, err) = (&mut std::io::stdout().lock(), &mut std::io::stderr().lock());
                ExitCode::from(match verb {
                    "outside-calls" => og_steps::pgx_rules::outside_calls_cli(&rest[1..], out, err),
                    "hla-concordance" => {
                        og_steps::pgx_rules::hla_concordance_cli(&rest[1..], out, err)
                    }
                    _ => og_steps::pypgx_summary::cli(&rest[1..], out, err),
                })
            }
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
        "dev" if rest.first().is_some_and(|v| v == "freshness") => {
            ExitCode::from(og_fresh::cli::main(&rest[1..]))
        }
        "dev" => ExitCode::from(og_dev::main(rest)),
        "recheck" => ExitCode::from(og_recheck::cli::main(rest)),
        "family" => ExitCode::from(og_steps::family::cli(rest)),
        "paraphase" => ExitCode::from(og_steps::paraphase::cli(
            rest,
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
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
        "traits" => ExitCode::from(og_steps::traits::cli(rest)),
        "bloodtype" => match rest.first() {
            Some(verb) => ExitCode::from(og_bloodtype::cli(
                verb,
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            None => usage_error("the following arguments are required: bloodtype tool"),
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
                "y" => og_reports::y::cli(verb, args, out, err),
                _ => og_reports::exomiser::cli(verb, args, out, err),
            })
        }
        "lookup" => ExitCode::from(og_lookup::cli(
            rest,
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
        "ancestry-fine" => ExitCode::from(og_ancestry::cli(
            rest,
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
        "archaic" if rest.first().is_some_and(|v| v == "compare") => {
            io(og_ancestry::compare::compare_cli, &rest[1..])
        }
        "archaic" if rest.first().is_some_and(|v| v == "genomes") => {
            io(og_ancestry::compare::genomes_cli, &rest[1..])
        }
        "archaic" => ExitCode::from(og_ancestry::archaic::cli(
            rest,
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
        "ancient" => ExitCode::from(og_ancestry::ancient::cli(
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
            Some("table") => ExitCode::from(og_summary::prs::cli(
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            _ => usage_error("og prs takes format or table"),
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
        "summary" => match rest.first().map(String::as_str) {
            Some("collect") => ExitCode::from(og_summary::collect::cli(
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            _ => usage_error("og summary takes collect"),
        },
        "report" => match rest.first().map(String::as_str) {
            Some("render") => io(
                |args, out, err| {
                    og_render::render::cli_with_spans(args, out, err, og_live::results::closable_of)
                },
                &rest[1..],
            ),
            Some("index") => io(
                |args, out, err| {
                    og_render::index::cli_using(args, out, err, og_live::results::closable_of)
                },
                &rest[1..],
            ),
            Some("charts") => io(og_render::charts::cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown report tool '{other}'")),
            None => usage_error("the following arguments are required: report tool"),
        },
        "live" => match rest.first().map(String::as_str) {
            Some("serve") => serve(og_live::server::cli, &rest[1..]),
            Some("results") => io(og_live::results::cli, &rest[1..]),
            Some("family") => io(og_live::family::cli, &rest[1..]),
            Some("demo") => io(og_live::demo::cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown live tool '{other}'")),
            None => usage_error("the following arguments are required: live tool"),
        },
        "view" => match rest.first().map(String::as_str) {
            Some("serve") => serve(og_live::viewer::cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown view tool '{other}'")),
            None => usage_error("the following arguments are required: view tool"),
        },
        "facts" => io(og_facts::cli, rest),
        "couple" => match rest.first().map(String::as_str) {
            Some("screen") => io(og_couple::screen::cli, &rest[1..]),
            Some("genes") => io(og_couple::genes::cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown couple tool '{other}'")),
            None => usage_error("the following arguments are required: couple tool"),
        },
        "longread" => match rest.first().map(String::as_str) {
            Some("clair3-model") => io(og_steps::longread::clair3_model_cli, &rest[1..]),
            Some(other) => usage_error(&format!("unknown longread tool '{other}'")),
            None => usage_error("the following arguments are required: longread tool"),
        },
        "roh" => match rest.first().map(String::as_str) {
            Some("summary") => ExitCode::from(og_summary::collect::group_c::roh::cli(
                &rest[1..],
                &mut std::io::stdout().lock(),
                &mut std::io::stderr().lock(),
            )),
            Some("af-file") => ExitCode::from(og_summary::collect::group_c::roh::af_file_cli(
                &rest[1..],
                &mut std::io::stderr().lock(),
            )),
            _ => usage_error("og roh takes summary or af-file"),
        },
        "mito" => ExitCode::from(og_summary::mito::cli(
            rest.first().map(String::as_str),
            rest.get(1..).unwrap_or_default(),
            &mut std::io::stdout().lock(),
            &mut std::io::stderr().lock(),
        )),
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
