//! `og pgx refblock-fill`: PharmCAT's multi-base positions called 0/0 from
//! the reference blocks of a gVCF (steps 07 and PHARMCAT_PREPROCESS run it).
//!
//! The fixture gVCF (`tests/fixtures/pharmcat/refblocks.g.vcf`) has
//! reference blocks over nine PharmCAT positions longer than one base, from
//! the subset of PharmCAT 3.4.0's positions file beside it. The test expands
//! the blocks into one-base records as `bcftools convert --gvcf2vcf |
//! bcftools view --trim-alt-alleles -i 'GT!="mis"'` does, runs the verb on
//! that stream in a child copy of this test binary ([`og_child`]), so its
//! stdin, stdout and stderr are the `og` binary's, and checks:
//!
//! 1. five positions are called 0/0, with PharmCAT's ID, REF and ALT, the
//!    flag PGX_REFBLOCK, the block's genotype (0 on haploid chrX) and the
//!    lowest GQ over the span: DPYD (two blocks meet inside the span),
//!    UGT1A1 *28, a TPMT insertion, a RYR1 two-base substitution and a G6PD
//!    deletion;
//! 2. four stay missing: a ./. hole inside the span (CYP2C9, 11 bp), a block
//!    with GQ 10, a variant record inside the span, a span no block reaches;
//! 3. each filled record takes the place of the one-base record at its
//!    position, every other line passes unchanged, and the output stays
//!    sorted;
//! 4. the GQ threshold: 20 by default, a block at exactly `--min-gq` counts;
//! 5. a block that stops one base short, a block with a failed FILTER, a ./.
//!    block with a high GQ, or a deletion that reaches into the span from
//!    before it, fills nothing; a gVCF with two samples or no #CHROM line is
//!    an error (exit 2).

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const CHILD: &str = "OG_PORT_TESTS_PGX_CHILD";
const MARK: &str = "\u{1e}og-child\u{1e}";

fn fx(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/fixtures/pharmcat")
        .join(name)
}

fn gvcf() -> String {
    fs::read_to_string(fx("refblocks.g.vcf")).unwrap()
}

/// (chrom, pos) -> (REF, ALT, GT, GQ, ID) of the five filled positions.
const FILLED: [((&str, u64), [&str; 5]); 5] = [
    (("chr1", 97_450_065), ["TG", "T", "0/0", "45", "rs72549303"]),
    (
        ("chr2", 233_760_233),
        ["CAT", "C,CATAT,CATATAT", "0/0", "50", "rs3064744"],
    ),
    (
        ("chr6", 18_149_032),
        ["C", "CT", "0/0", "40", "rs759836180"],
    ),
    (
        ("chr19", 38_572_266),
        ["TC", "CT", "0/0", "35", "rs193922862"],
    ),
    (("chrX", 154_532_055), ["CTCT", "C", "0", "30", "."]),
];
const LEFT: [((&str, u64), &str); 4] = [
    (("chr10", 94_942_205), "a ./. hole inside the REF span"),
    (("chr10", 94_949_161), "a block with GQ 10"),
    (
        ("chr10", 94_949_281),
        "a variant record inside the REF span",
    ),
    (("chr6", 18_149_103), "no block over the span"),
];
const UGT1A1: (&str, u64) = ("chr2", 233_760_233);

/// Not a test of its own: the child [`run`] starts runs `og pgx
/// refblock-fill` here with the argv it was handed, then exits with the
/// verb's code. In the normal test run it returns at once.
#[test]
fn og_child() {
    let Some(argv) = std::env::var_os(CHILD) else {
        return;
    };
    let argv: Vec<String> = argv
        .to_string_lossy()
        .split('\u{1f}')
        .map(str::to_owned)
        .collect();
    print!("{MARK}");
    eprint!("{MARK}");
    let _ = std::io::stdout().flush();
    let code = og_filters::refblock::cli(&argv);
    let _ = std::io::stdout().flush();
    let _ = std::io::stderr().flush();
    std::process::exit(i32::from(code));
}

/// What gvcf2vcf plus `view --trim-alt-alleles -i 'GT!="mis"'` writes: one
/// record per base of each block (REF from the reference: N here, the verb
/// reads only its length), ALT . for a reference call, no ./. records.
fn expand(gvcf: &str) -> String {
    let mut out = String::new();
    for line in gvcf.split_inclusive('\n') {
        if line.starts_with('#') {
            out.push_str(line);
            continue;
        }
        let c: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        let gt = c[9].split(':').next().unwrap();
        if gt.split(['/', '|']).any(|a| a == ".") {
            continue;
        }
        let alts: Vec<&str> = c[4].split(',').filter(|a| *a != "<*>").collect();
        let end = c[7].find("END=").and_then(|i| {
            let digits: String = c[7][i + 4..]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect();
            digits.parse::<u64>().ok()
        });
        match end {
            Some(end) if alts.is_empty() => {
                for p in c[1].parse::<u64>().unwrap()..=end {
                    let p = p.to_string();
                    out.push_str(
                        &[c[0], &p, ".", "N", ".", c[5], c[6], ".", c[8], c[9]].join("\t"),
                    );
                    out.push('\n');
                }
            }
            _ => {
                let alt = if alts.is_empty() {
                    ".".to_owned()
                } else {
                    alts.join(",")
                };
                let mut cols: Vec<&str> = c[..4].to_vec();
                cols.push(&alt);
                cols.extend(&c[5..]);
                out.push_str(&cols.join("\t"));
                out.push('\n');
            }
        }
    }
    out
}

static RUNS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// The verb on `gvcf` (its expansion on stdin, or `stdin`): exit code,
/// stdout, stderr.
fn run(gvcf: &str, opts: &[&str], stdin: Option<&str>) -> (i32, String, String) {
    let n = RUNS.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("pgx-refblock-{}-{n}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let g = dir.join("records.vcf");
    fs::write(&g, gvcf).unwrap();
    let input = dir.join("stdin.vcf");
    fs::write(&input, stdin.map_or_else(|| expand(gvcf), str::to_owned)).unwrap();
    let mut argv: Vec<String> = opts.iter().map(|s| (*s).to_owned()).collect();
    argv.push(
        fx("pharmcat_positions-3.4.0-subset.vcf")
            .to_string_lossy()
            .into_owned(),
    );
    argv.push(g.to_string_lossy().into_owned());
    let out = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "og_child",
            "--nocapture",
            "--test-threads=1",
            "-q",
        ])
        .env(CHILD, argv.join("\u{1f}"))
        .stdin(Stdio::from(File::open(&input).unwrap()))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap()
        .wait_with_output()
        .unwrap();
    let _ = fs::remove_dir_all(&dir);
    let after = |b: &[u8]| {
        let s = String::from_utf8_lossy(b);
        match s.split_once(MARK) {
            Some((_, rest)) => rest.to_owned(),
            None => panic!("the child never reached the verb: {s}"),
        }
    };
    (
        out.status.code().unwrap_or(-1),
        after(&out.stdout),
        after(&out.stderr),
    )
}

/// Stdout of a run on valid input, which must exit 0, so "nothing filled"
/// never comes from a run that failed.
fn ok(gvcf: &str, opts: &[&str]) -> String {
    let (rc, out, err) = run(gvcf, opts, None);
    assert_eq!(rc, 0, "a valid input exits 0: {err}");
    out
}

fn records(vcf: &str) -> Vec<Vec<String>> {
    vcf.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

fn key(r: &[String]) -> (String, u64) {
    (r[0].clone(), r[1].parse().unwrap())
}

fn fills(vcf: &str) -> BTreeMap<(String, u64), Vec<String>> {
    records(vcf)
        .into_iter()
        .filter(|r| r[7] == "PGX_REFBLOCK")
        .map(|r| (key(&r), r))
        .collect()
}

fn filled(vcf: &str, (chrom, pos): (&str, u64)) -> bool {
    fills(vcf).contains_key(&(chrom.to_owned(), pos))
}

#[test]
fn the_fixture_fills_five_and_leaves_four() {
    let g = gvcf();
    let (rc, out, err) = run(&g, &[], None);
    assert_eq!(rc, 0, "runs and exits 0: {err}");
    let got = fills(&out);
    // 1. filled
    let mut want: Vec<(String, u64)> = FILLED
        .iter()
        .map(|((c, p), _)| ((*c).to_owned(), *p))
        .collect();
    want.sort();
    assert_eq!(
        got.keys().cloned().collect::<Vec<_>>(),
        want,
        "five positions called 0/0: DPYD, UGT1A1 *28, TPMT, RYR1, G6PD"
    );
    for ((chrom, pos), [r, alt, gt, gq, id]) in FILLED {
        let row = got.get(&(chrom.to_owned(), pos));
        assert!(
            row.is_some_and(|x| x[2..5] == [id, r, alt]
                && x[6] == "PASS"
                && x[8] == "GT:GQ"
                && x[9] == format!("{gt}:{gq}")),
            "{chrom}:{pos}: PharmCAT's ID, REF and ALT, GT {gt}, lowest GQ {gq}: {row:?}"
        );
    }
    // 2. left missing
    for ((chrom, pos), why) in LEFT {
        assert!(
            !got.contains_key(&(chrom.to_owned(), pos)),
            "{chrom}:{pos} stays missing ({why})"
        );
    }
    assert!(
        err.contains(
            "5 of 9 PharmCAT positions longer than one base called 0/0 from reference blocks with GQ >= 20"
        ),
        "the summary on stderr counts 5 of the 9 positions longer than one base: {err}"
    );
    // 3. replaced in place, the rest unchanged, sorted
    let inp = records(&expand(&g));
    let outr = records(&out);
    let mut at: BTreeMap<(String, u64), usize> = BTreeMap::new();
    for r in &outr {
        *at.entry(key(r)).or_default() += 1;
    }
    let filled_keys: Vec<(String, u64)> = FILLED
        .iter()
        .map(|((c, p), _)| ((*c).to_owned(), *p))
        .collect();
    assert!(
        filled_keys.iter().all(|k| at.get(k) == Some(&1)),
        "one record at each filled position (the one-base record is replaced)"
    );
    let kept: Vec<&Vec<String>> = inp
        .iter()
        .filter(|r| !filled_keys.contains(&key(r)))
        .collect();
    let passed: Vec<&Vec<String>> = outr.iter().filter(|r| r[7] != "PGX_REFBLOCK").collect();
    assert_eq!(
        passed, kept,
        "every other record passes unchanged and in order"
    );
    assert!(
        at.contains_key(&("chr2".to_owned(), 233_760_234))
            && at.contains_key(&("chr2".to_owned(), 233_760_235)),
        "the bases after a filled position inside its span stay (chr2:233760234-233760235)"
    );
    let order = ["chr1", "chr2", "chr6", "chr10", "chr19", "chrX"];
    let keys: Vec<(usize, u64)> = outr
        .iter()
        .map(|r| {
            (
                order.iter().position(|c| *c == r[0]).unwrap(),
                r[1].parse().unwrap(),
            )
        })
        .collect();
    let mut sorted = keys.clone();
    sorted.sort_unstable();
    assert_eq!(keys, sorted, "the output is sorted");
    let head: Vec<&str> = out.lines().filter(|l| l.starts_with('#')).collect();
    assert!(
        head.len() >= 2
            && head[head.len() - 1].starts_with("#CHROM")
            && head[head.len() - 2].starts_with("##INFO=<ID=PGX_REFBLOCK,"),
        "the header declares PGX_REFBLOCK before #CHROM: {head:?}"
    );
    assert_eq!(
        head.last()
            .map(|h| h.split('\t').skip(9).collect::<Vec<_>>()),
        Some(vec!["REFBLOCKS"]),
        "the sample column stays"
    );
}

#[test]
fn the_gq_threshold() {
    let g = gvcf();
    let o = ok(&g, &["--min-gq", "50"]);
    assert!(
        filled(&o, UGT1A1),
        "--min-gq 50: UGT1A1 (GQ 50) still filled, a block at exactly the threshold counts"
    );
    assert!(
        !filled(&o, ("chr1", 97_450_065)),
        "--min-gq 50: DPYD (one block GQ 45) no longer filled"
    );
    let o = ok(&g, &["--min-gq", "51"]);
    assert!(fills(&o).is_empty(), "--min-gq 51: nothing filled");
    // the first GQ-50 block of a filled position: UGT1A1
    let low = g.replacen("0/0:50:30\n", "0/0:19:30\n", 1);
    let o = ok(&low, &[]);
    assert!(
        !filled(&o, UGT1A1),
        "default threshold 20: UGT1A1 with a GQ 19 block is not filled"
    );
}

#[test]
fn edges_that_fill_nothing_and_inputs_that_are_refused() {
    let g = gvcf();
    let hole = g.replace("./.:0:0", "./.:50:0");
    let o = ok(&hole, &[]);
    assert!(
        hole != g && !filled(&o, ("chr10", 94_942_205)),
        "the ./. hole stays a hole when its GQ is 50 (the genotype decides, not only the GQ)"
    );
    // A deletion that starts before the span and reaches into it, with a
    // block over the deleted bases after it (as GATK writes them): the
    // block covers the span, the deletion still overlaps it.
    let dele = g.replace(
        "chr2\t233759900\t.\tT\t<*>\t0\t.\tEND=233761000\tGT:GQ:MIN_DP\t0/0:50:30\n",
        "chr2\t233759900\t.\tT\t<*>\t0\t.\tEND=233760230\tGT:GQ:MIN_DP\t0/0:50:30\n\
         chr2\t233760231\t.\tGCCA\tG,<*>\t30\tPASS\t.\tGT:GQ:DP\t0/1:30:28\n\
         chr2\t233760232\t.\tC\t<*>\t0\t.\tEND=233761000\tGT:GQ:MIN_DP\t0/0:50:30\n",
    );
    let o = ok(&dele, &[]);
    assert!(
        dele != g && !filled(&o, UGT1A1),
        "a deletion from before the span that overlaps it fills nothing, though blocks cover the span"
    );
    let short = g.replace("END=233761000", "END=233760234");
    let o = ok(&short, &[]);
    assert!(
        !filled(&o, UGT1A1),
        "a block ending one base before the span ends fills nothing"
    );
    let lowq = g.replace(
        "T\t<*>\t0\t.\tEND=233761000",
        "T\t<*>\t0\tLowQual\tEND=233761000",
    );
    let o = ok(&lowq, &[]);
    assert!(
        lowq != g && !filled(&o, UGT1A1),
        "a block with a failed FILTER fills nothing"
    );
    let two = g.replace("\tREFBLOCKS\n", "\tREFBLOCKS\tOTHER\n");
    let (rc, _, e) = run(&two, &[], None);
    assert!(
        rc == 2 && e.contains("one sample"),
        "a gVCF with two samples: exit 2: {rc} {e}"
    );
    let nohead: String = g
        .split_inclusive('\n')
        .filter(|l| !l.starts_with("#CHROM"))
        .collect();
    let (rc, _, e) = run(&nohead, &[], Some(&expand(&g)));
    assert!(
        rc == 2 && e.contains("#CHROM"),
        "a gVCF without its #CHROM line: exit 2: {rc} {e}"
    );
}
