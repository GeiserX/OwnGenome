//! `og dev demo-sample`: DEMO-001, the made-up sample of the docs pictures
//! and the live-site demo.
//!
//! ```text
//! og dev demo-sample --genome-dir DIR [--seed 1] [--records 4700000]
//!                    [--sample DEMO-001] [--sex male|female] [--repo DIR] [--og OG]
//! og dev demo-sample check --genome-dir DIR [--sample DEMO-001]
//! ```
//!
//! The first form writes `DIR/<sample>` with every file steps 24, 27 and 28
//! read, as the steps that make them would have written them, so every card
//! of the reports shows a value. Nothing is read from a real sample: the
//! values come from a small seeded generator ([`Rng`]) and from a few
//! hand-picked numbers. It also writes the three reference files the run
//! manifest reads (`DIR/clinvar`, `DIR/t1k_idx`, `DIR/prs_scores`), cut down
//! to their header lines.
//!
//! Every table a pipeline step writes with og is written here by the same og
//! verb (`--og`, default this binary), fed made-up inputs, so the demo always
//! shows what the product writes: `og pgx outside-calls` and
//! `og pgx cpic-report` (steps 36 and 27), `og prs table` (25),
//! `og qc sample` and `og qc sex-loss` (33), `og ancestry-fine report`
//! (26b), `og paraphase table` (35b), `og couple screen` (39),
//! `og traits table` (42) and `og ancient report` (44). The other tables are
//! written in their writers' columns, with the values those writers give for
//! the made-up inputs named in each function. Each part lives in the file of
//! its report chapter: [`qc`], [`clinical`], [`clinical_clinvar`], [`pgx`],
//! [`sv`], [`lineage`] and [`traits`].
//!
//! The generator is not Python's: the replaced `make_demo_sample.py` drew
//! from Python's Mersenne Twister, and the values here differ from its
//! values. What stays is every hand-picked number and every planted finding
//! (the genes and positions of [`clinical_clinvar::PLANTED`] and
//! [`FINDINGS`]), so every check that names them still holds.
//!
//! `--sample` and `--sex` write another invented person the same way: the
//! live-site demo writes DEMO-002, DEMO-001's partner (female), and
//! DEMO-003, their son. A female gets no Y files and her sex checks read two
//! X. DEMO-002 carries an FMR1 allele in Stranger's premutation range, a
//! CHIP call and a mosaic loss; DEMO-003 an AnnotSV class 4 deletion and the
//! planted results to act on.
//!
//! The second form reads `DIR/<sample>/summary.json` (step 24 writes it) and
//! fails when a section is not "ok" or the CYP2D6 callers disagree, so a new
//! report card the demo does not feed yet is caught before a picture shows
//! "Not run".
//!
//! Exit 0 when written (or every section is ok), 1 when a part fails or a
//! section is not ok, 2 on a usage error.

pub mod ancestry_paint;
pub mod clinical;
pub mod clinical_clinvar;
pub mod family;
pub mod lineage;
pub mod pgx;
pub mod qc;
pub mod sv;
pub mod traits;

use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const USAGE: &str = "usage: og dev demo-sample --genome-dir DIR [--seed N] [--records N] [--sample DEMO-001] [--sex male|female] [--repo DIR] [--og OG]\n       og dev demo-sample check --genome-dir DIR [--sample DEMO-001]\n";

/// The sample the docs pictures show.
pub const SAMPLE: &str = "DEMO-001";

/// GRCh38 primary chromosomes and their lengths.
pub const CHROMS: [(&str, i64); 24] = [
    ("chr1", 248956422),
    ("chr2", 242193529),
    ("chr3", 198295559),
    ("chr4", 190214555),
    ("chr5", 181538259),
    ("chr6", 170805979),
    ("chr7", 159345973),
    ("chr8", 145138636),
    ("chr9", 138394717),
    ("chr10", 133797422),
    ("chr11", 135086622),
    ("chr12", 133275309),
    ("chr13", 114364328),
    ("chr14", 107043718),
    ("chr15", 101991189),
    ("chr16", 90338345),
    ("chr17", 83257441),
    ("chr18", 80373285),
    ("chr19", 58617616),
    ("chr20", 64444167),
    ("chr21", 46709983),
    ("chr22", 50818468),
    ("chrX", 156040895),
    ("chrY", 57227415),
];

/// chr1 to chr22.
pub fn autosomes() -> &'static [(&'static str, i64)] {
    &CHROMS[..22]
}

/// Where `chrom` sits in [`CHROMS`] (unknown names last).
pub fn chrom_index(chrom: &str) -> usize {
    CHROMS
        .iter()
        .position(|(c, _)| *c == chrom)
        .unwrap_or(CHROMS.len())
}

/// The mean read depth of the invented genome.
pub const MEAN_DEPTH: f64 = 31.3;

/// When the invented run started: 2026-10-01 09:00 UTC, the time of the
/// demo's PharmCAT report too (`logs/run_status.tsv`, the report's "Latest
/// run" line). Fixed, so every retake shows the same line; the files written
/// now are newer, so none is marked stale.
pub const RUN_STARTED: i64 = 1_790_845_200;

/// The invented findings of the demo's other people, one of each live-site
/// widget's marks.
pub const FINDINGS: &[(&str, &[&str])] = &[
    ("DEMO-002", &["fmr1", "chip", "mocha"]),
    ("DEMO-003", &["annotsv", "act"]),
];

/// What one run of the writer knows: who, where, and with which og.
pub struct Ctx {
    pub repo: PathBuf,
    pub og: PathBuf,
    pub genome_dir: PathBuf,
    /// `genome_dir/<sample>`.
    pub dir: PathBuf,
    pub sample: String,
    pub male: bool,
    pub records: usize,
}

impl Ctx {
    /// "male" or "female".
    pub fn sex(&self) -> &'static str {
        if self.male { "male" } else { "female" }
    }

    /// Does this person carry the invented finding `what` ([`FINDINGS`])?
    pub fn has(&self, what: &str) -> bool {
        FINDINGS
            .iter()
            .any(|(s, f)| *s == self.sample && f.contains(&what))
    }

    /// The chromosomes the person has reads and variants on: a female has no chrY.
    pub fn own_chroms(&self) -> Vec<(&'static str, i64)> {
        CHROMS
            .iter()
            .copied()
            .filter(|(c, _)| self.male || *c != "chrY")
            .collect()
    }

    /// Is `chrom` a single copy here (X and Y of a male)?
    pub fn haploid(&self, chrom: &str) -> bool {
        self.male && (chrom == "chrX" || chrom == "chrY")
    }

    /// `dir/<rel>`.
    pub fn path(&self, rel: &str) -> PathBuf {
        self.dir.join(rel)
    }

    /// A file of the repository (`assets/...`).
    pub fn asset(&self, rel: &str) -> PathBuf {
        self.repo.join(rel)
    }

    /// A VCF header over the GRCh38 contigs, with `extra` lines and the
    /// sample column.
    pub fn vcf_header(&self, extra: &[&str]) -> String {
        let mut s = String::from("##fileformat=VCFv4.2\n");
        for (c, n) in CHROMS {
            s.push_str(&format!("##contig=<ID={c},length={n}>\n"));
        }
        for x in extra {
            s.push_str(x);
            s.push('\n');
        }
        s.push_str(&format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            self.sample
        ));
        s
    }

    /// Runs `og ARGS`, its stdout dropped; an error names the verb and its stderr.
    pub fn og(&self, args: &[&str]) -> Result<(), String> {
        self.og_output(args).map(|_| ())
    }

    /// Runs `og ARGS`; its stdout, or an error naming the verb and its stderr.
    pub fn og_output(&self, args: &[&str]) -> Result<String, String> {
        let verb = args.iter().take(2).copied().collect::<Vec<_>>().join(" ");
        let out = Command::new(&self.og)
            .args(args)
            .stdin(Stdio::null())
            .output()
            .map_err(|e| format!("cannot run {}: {e}", self.og.display()))?;
        if !out.status.success() {
            let code = out
                .status
                .code()
                .map_or_else(|| "a signal".to_owned(), |c| c.to_string());
            return Err(format!(
                "og {verb} failed (exit {code}): {}",
                String::from_utf8_lossy(&out.stderr).trim()
            ));
        }
        String::from_utf8(out.stdout).map_err(|_| format!("og {verb}: stdout is not UTF-8"))
    }
}

/// A path as an argument.
pub fn arg(p: &Path) -> &str {
    p.to_str().unwrap_or_default()
}

/// Writes `text` to `path`, making its folder.
pub fn put(path: &Path, text: &str) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// Writes `text` gzipped to `path` (level 6, no time stamp, so a retake
/// writes the same bytes).
pub fn put_gz(path: &Path, text: &str) -> Result<(), String> {
    let mut gz = gz_writer(path, 6)?;
    gz.write_all(text.as_bytes())
        .map_err(|e| format!("{}: {e}", path.display()))?;
    finish_gz(gz, path)
}

/// A gzip writer on a new file at `path`.
pub fn gz_writer(
    path: &Path,
    level: u32,
) -> Result<GzEncoder<std::io::BufWriter<fs::File>>, String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    let f = fs::File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(flate2::GzBuilder::new()
        .mtime(0)
        .write(std::io::BufWriter::new(f), Compression::new(level)))
}

/// Closes a [`gz_writer`].
pub fn finish_gz(gz: GzEncoder<std::io::BufWriter<fs::File>>, path: &Path) -> Result<(), String> {
    gz.finish()
        .and_then(|mut w| w.flush())
        .map_err(|e| format!("{}: {e}", path.display()))
}

/// Reads a text file.
pub fn read(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

/// Copies `from` to `to`, making the folder of `to`.
pub fn copy(from: &Path, to: &Path) -> Result<(), String> {
    if let Some(parent) = to.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    fs::copy(from, to)
        .map(|_| ())
        .map_err(|e| format!("{} -> {}: {e}", from.display(), to.display()))
}

/// The rows of a tab-separated table after its `#` comment lines and blank
/// lines: the header's names zipped with each row (a short row gets "").
pub fn table(text: &str) -> Vec<Vec<(String, String)>> {
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let Some(head) = lines.next() else {
        return Vec::new();
    };
    let names: Vec<&str> = head.split('\t').collect();
    lines
        .map(|l| {
            let cells: Vec<&str> = l.split('\t').collect();
            names
                .iter()
                .enumerate()
                .map(|(i, n)| {
                    (
                        (*n).to_owned(),
                        cells.get(i).copied().unwrap_or("").to_owned(),
                    )
                })
                .collect()
        })
        .collect()
}

/// A cell of a [`table`] row ("" when the column is missing).
pub fn cell<'a>(row: &'a [(String, String)], name: &str) -> &'a str {
    row.iter()
        .find(|(n, _)| n == name)
        .map_or("", |(_, v)| v.as_str())
}

/// The value of `VAR=` in `versions.env`, before any `#` comment, unquoted.
pub fn versions_env(repo: &Path, var: &str) -> Result<String, String> {
    let text = read(&repo.join("versions.env"))?;
    text.lines()
        .find_map(|l| l.strip_prefix(&format!("{var}=")))
        .map(|v| {
            v.split('#')
                .next()
                .unwrap_or("")
                .trim()
                .trim_matches('"')
                .to_owned()
        })
        .ok_or_else(|| format!("{var} is not in versions.env"))
}

/// A scratch folder under the system's temporary folder, removed when dropped.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new(name: &str) -> Result<Scratch, String> {
        let base = std::env::temp_dir();
        for k in 0..1000u32 {
            let p = base.join(format!("og-demo-{}-{name}-{k}", std::process::id()));
            if fs::create_dir(&p).is_ok() {
                return Ok(Scratch(p));
            }
        }
        Err(format!(
            "cannot make a scratch folder in {}",
            base.display()
        ))
    }

    pub fn join(&self, rel: &str) -> PathBuf {
        self.0.join(rel)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// `epoch` as `YYYY-MM-DD HH:MM:SS UTC`.
pub fn utc(epoch: i64) -> String {
    let days = epoch.div_euclid(86_400);
    let secs = epoch.rem_euclid(86_400);
    // Howard Hinnant's civil_from_days.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02} {:02}:{:02}:{:02} UTC",
        secs / 3600,
        secs % 3600 / 60,
        secs % 60
    )
}

/// The complementary error function, to about 1.2e-7 relative error
/// (Numerical Recipes' `erfcc`): enough for the two or three decimals the
/// coverage curves are written with.
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 { r } else { 2.0 - r }
}

/// The demo's random numbers: xoshiro256** seeded through SplitMix64 from a
/// text seed. Each part gets its own seed (`"<seed>:<part>"`), so changing
/// one part leaves the others' values alone. Not Python's generator: the
/// demo is made up, so nothing is gained by drawing Python's numbers.
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn new(seed: &str) -> Rng {
        // FNV-1a of the seed text, spread over the state by SplitMix64.
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for b in seed.bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
        let mut sm = h;
        let mut next = || {
            sm = sm.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = sm;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        };
        Rng {
            s: [next(), next(), next(), next()],
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// A float in [0, 1).
    pub fn random(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// An integer in [0, n), with no bias.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "below(0)");
        let zone = u64::MAX - u64::MAX % n;
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % n;
            }
        }
    }

    /// An integer in [lo, hi).
    pub fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + self.below((hi - lo) as u64) as i64
    }

    /// An integer in [lo, hi], both ends included.
    pub fn int(&mut self, lo: i64, hi: i64) -> i64 {
        self.range(lo, hi + 1)
    }

    /// An index into a slice of `n`.
    pub fn index(&mut self, n: usize) -> usize {
        self.below(n as u64) as usize
    }

    /// A float in [a, b).
    pub fn uniform(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.random()
    }

    /// A normal draw (Box-Muller).
    pub fn gauss(&mut self, mu: f64, sigma: f64) -> f64 {
        let u1 = 1.0 - self.random();
        let u2 = self.random();
        mu + sigma * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// One item of `xs`.
    pub fn pick<'a, T>(&mut self, xs: &'a [T]) -> &'a T {
        &xs[self.index(xs.len())]
    }

    /// An index drawn with the given weights.
    pub fn weighted(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut x = self.random() * total;
        for (i, w) in weights.iter().enumerate() {
            if x < *w {
                return i;
            }
            x -= w;
        }
        weights.len() - 1
    }

    /// `k` different indices of `0..n`, in draw order.
    pub fn sample_idx(&mut self, n: usize, k: usize) -> Vec<usize> {
        assert!(k <= n, "sample of {k} from {n}");
        let mut pool: Vec<usize> = (0..n).collect();
        for i in 0..k {
            let j = i + self.index(n - i);
            pool.swap(i, j);
        }
        pool.truncate(k);
        pool
    }

    /// Two different letters of `ACGT`.
    pub fn two_bases(&mut self) -> (char, char) {
        let i = self.sample_idx(4, 2);
        (BASES[i[0]] as char, BASES[i[1]] as char)
    }

    /// `k` different integers of [lo, hi), sorted.
    pub fn distinct_sorted(&mut self, lo: i64, hi: i64, k: usize) -> Vec<i64> {
        assert!((hi - lo) as usize >= k, "{k} positions in [{lo}, {hi})");
        let mut out: Vec<i64> = Vec::with_capacity(k);
        while out.len() < k {
            let missing = k - out.len();
            out.extend((0..missing).map(|_| self.range(lo, hi)));
            out.sort_unstable();
            out.dedup();
        }
        out
    }
}

/// The four bases.
pub const BASES: &[u8; 4] = b"ACGT";

/// One part of the sample: its name (and seed) and its writer.
type Part = (&'static str, fn(&Ctx, &mut Rng) -> Result<(), String>);

/// The parts in the order they are written: `read_count` and `mosaic` read
/// the files of `reads` and `sample_qc`, so they come after them.
const PARTS: &[Part] = &[
    ("run_status", qc::run_status),
    ("variants", qc::variants),
    ("clinvar", clinical_clinvar::clinvar),
    ("pgx", pgx::pgx),
    ("prs", lineage::prs),
    ("coverage", qc::coverage),
    ("sample_qc", qc::sample_qc),
    ("reads", qc::reads),
    ("structural", sv::structural),
    ("expansions", clinical::expansions),
    ("telomere", traits::telomere),
    ("roh", lineage::roh),
    ("mito", clinical::mito),
    ("haplogroups", lineage::haplogroups),
    ("clinical", clinical::clinical),
    ("read_count", qc::read_count),
    ("mosaic", qc::mosaic),
    ("ancestry_fine", lineage::ancestry_fine),
    ("paralogs", sv::paralogs),
    ("couple", clinical_clinvar::couple),
    ("consumer_export", traits::consumer_export),
    ("bloodtype", traits::bloodtype),
    ("traits", traits::traits),
    ("y_str", lineage::y_str),
    ("viral", traits::viral),
    ("lpa", clinical::lpa),
    ("chip", clinical::chip),
    ("archaic", lineage::archaic),
    ("ancient", lineage::ancient),
    ("exomiser", clinical::exomiser),
    ("mocha", clinical::mocha),
    ("vntyper", clinical::vntyper),
    ("mei", sv::mei),
    ("locityper", pgx::locityper),
    ("callability", qc::callability),
];

/// The repository holding `assets/` and `versions.env`: the current folder
/// or one above it.
fn find_repo() -> Option<PathBuf> {
    let mut dir = std::env::current_dir().ok()?;
    loop {
        if dir.join("assets").join("steps.tsv").is_file() && dir.join("versions.env").is_file() {
            return Some(dir);
        }
        if !dir.pop() {
            return None;
        }
    }
}

/// `og dev demo-sample ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let usage = |err: &mut dyn Write, msg: &str| {
        let _ = writeln!(err, "{USAGE}og dev demo-sample: error: {msg}");
        2
    };
    let (checking, rest) = match args.first().map(String::as_str) {
        Some("check") => (true, &args[1..]),
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        _ => (false, args),
    };
    let mut opts: Vec<(&str, &str)> = Vec::new();
    let mut i = 0;
    while i < rest.len() {
        let name = rest[i].as_str();
        let allowed: &[&str] = if checking {
            &["--genome-dir", "--sample"]
        } else {
            &[
                "--genome-dir",
                "--seed",
                "--records",
                "--sample",
                "--sex",
                "--repo",
                "--og",
            ]
        };
        if !allowed.contains(&name) {
            return usage(err, &format!("unrecognized argument: {name}"));
        }
        let Some(value) = rest.get(i + 1) else {
            return usage(err, &format!("{name} needs a value"));
        };
        opts.push((name, value.as_str()));
        i += 2;
    }
    let get = |n: &str| opts.iter().rev().find(|(k, _)| *k == n).map(|(_, v)| *v);
    let Some(genome_dir) = get("--genome-dir") else {
        return usage(err, "--genome-dir is required");
    };
    let sample = get("--sample").unwrap_or(SAMPLE);
    if checking {
        return check(Path::new(genome_dir), sample, out, err);
    }
    let valid_name = sample.len() == 8
        && sample.starts_with("DEMO-")
        && sample[5..].bytes().all(|b| b.is_ascii_digit());
    if !valid_name {
        return usage(
            err,
            &format!("--sample {sample}: an invented name of the form DEMO-001"),
        );
    }
    let male = match get("--sex").unwrap_or("male") {
        "male" => true,
        "female" => false,
        other => {
            return usage(
                err,
                &format!("--sex: invalid choice: '{other}' (choose from 'male', 'female')"),
            );
        }
    };
    let seed = get("--seed").unwrap_or("1");
    if seed.parse::<i64>().is_err() {
        return usage(err, &format!("--seed: invalid int value: '{seed}'"));
    }
    let records = match get("--records").unwrap_or("4700000").parse::<usize>() {
        Ok(n) => n,
        Err(_) => return usage(err, "--records: a whole number of records"),
    };
    let repo = match get("--repo").map(PathBuf::from).or_else(find_repo) {
        Some(r) => r,
        None => return usage(err, "no assets/steps.tsv here or above: give --repo"),
    };
    let og = match get("--og") {
        Some(p) => PathBuf::from(p),
        None => match std::env::current_exe() {
            Ok(p) => p,
            Err(e) => return usage(err, &format!("cannot find this og: {e}; give --og")),
        },
    };
    let genome_dir = PathBuf::from(genome_dir);
    let dir = genome_dir.join(sample);
    if dir.exists() {
        let _ = writeln!(
            err,
            "og dev demo-sample: {} exists: give an empty --genome-dir",
            dir.display()
        );
        return 1;
    }
    if let Err(e) = fs::create_dir_all(&dir) {
        let _ = writeln!(err, "og dev demo-sample: {}: {e}", dir.display());
        return 1;
    }
    let ctx = Ctx {
        repo,
        og,
        genome_dir,
        dir,
        sample: sample.to_owned(),
        male,
        records,
    };
    for (name, part) in PARTS {
        let mut rng = Rng::new(&format!("{seed}:{name}"));
        if let Err(e) = part(&ctx, &mut rng) {
            let _ = writeln!(err, "og dev demo-sample: {name}: {e}");
            return 1;
        }
        let _ = writeln!(out, "  {name}: written");
    }
    let _ = writeln!(out, "{} written in {}", ctx.sample, ctx.dir.display());
    0
}

/// `og dev demo-sample check`: fail when a section of step 24's summary is
/// not ok, or the CYP2D6 callers do not agree.
fn check(genome_dir: &Path, sample: &str, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    use pycompat::json::Value;
    let path = genome_dir.join(sample).join("summary.json");
    let doc = match read(&path)
        .and_then(|t| pycompat::json::loads(&t).map_err(|e| format!("{}: {e}", path.display())))
    {
        Ok(d) => d,
        Err(e) => {
            let _ = writeln!(err, "og dev demo-sample check: {e}");
            return 1;
        }
    };
    let field = |v: &Value, k: &str| -> Option<Value> {
        match v {
            Value::Object(m) => m.get(k).cloned(),
            _ => None,
        }
    };
    let text = |v: Option<Value>| match v {
        Some(Value::Str(s)) => s,
        _ => String::new(),
    };
    let Some(Value::Object(sections)) = field(&doc, "sections") else {
        let _ = writeln!(
            err,
            "og dev demo-sample check: {} has no sections",
            path.display()
        );
        return 1;
    };
    let mut bad: Vec<(String, String, String)> = Vec::new();
    for (k, v) in &sections {
        let state = text(field(v, "state"));
        let _ = writeln!(out, "  {k:<14} {state}");
        if state != "ok" {
            bad.push((k.clone(), state, text(field(v, "note"))));
        }
    }
    let cyp = field(&doc, "cyp2d6").unwrap_or(Value::Null);
    if field(&cyp, "agree") != Some(Value::Bool(true)) {
        let calls = field(&cyp, "calls").unwrap_or(Value::Null);
        let calls =
            pycompat::json::dumps(&calls, &pycompat::json::Dump::default()).unwrap_or_default();
        bad.push((
            "cyp2d6".to_owned(),
            "callers do not agree".to_owned(),
            calls,
        ));
    }
    if !bad.is_empty() {
        for (k, st, note) in bad {
            let _ = writeln!(err, "{}", format!("FAIL: {k} is {st} {note}").trim_end());
        }
        let _ = writeln!(
            err,
            "A card would say 'Not run' or worse: feed that section in rust/crates/og-dev/src/demo/."
        );
        return 1;
    }
    let _ = writeln!(out, "Every section of {} is ok.", path.display());
    0
}
