//! `og dev scan images` and `og dev versions coupled`: the two checks of
//! `scripts/ci/check-images.sh` that were Python. The script keeps the
//! lists (which paths are scanned, which files and single literals are
//! exempt, with the reason) and passes them in.
//!
//! ```text
//! og dev scan images --root DIR --exempt-files TEXT --exempt-pairs TEXT < FILES
//! og dev versions coupled --env versions.env
//! ```
//!
//! `scan images` reads the files to scan (one path relative to DIR per
//! line) on stdin and prints one `FAIL: <file>:<line> names the image
//! <image>; ...` per literal image, then `scanned N files`. An image is
//! `registry/name:tag`, `name@sha256:...`, or a Docker Hub official name
//! (those `versions.env` uses and the usual base images) with a tag. A line
//! of `--exempt-pairs` (`file image reason`) allows one occurrence of that
//! image in that file; listing a pair twice allows two.
//!
//! `versions coupled` checks that `MINIMAP2_SAMTOOLS_IMAGE`, the one image
//! of the Nextflow alignment task, names in its `versions.env` comment the
//! minimap2 and samtools releases of `MINIMAP2_IMAGE` and `SAMTOOLS_IMAGE`.
//!
//! Exit 0 clean, 1 on a finding, 2 on a usage error or a `versions.env`
//! the strict reader refuses.
//!
//! Where it differs from the Python on purpose: an exempt pair whose file
//! no longer names its image as often as the list allows fails
//! (`FAIL: <file> no longer names <image>; drop its line ...`), so an
//! exemption cannot outlive its reason; and `versions.env` is read strictly
//! ([`crate::versions_env`]).

use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::LazyLock;

use regex::Regex;

use crate::lookaround::Rule;
use crate::versions_env::{self, Pin};

const USAGE: &str = "usage: og dev scan images --root DIR --exempt-files TEXT --exempt-pairs TEXT < FILES\n       og dev versions coupled --env versions.env\n";

/// Docker Hub official images that `versions.env` may not use today.
const BASES: [&str; 15] = [
    "ubuntu",
    "debian",
    "alpine",
    "busybox",
    "centos",
    "rockylinux",
    "fedora",
    "python",
    "node",
    "perl",
    "r-base",
    "openjdk",
    "eclipse-temurin",
    "golang",
    "rust",
];

/// Any tag: a version, or a word such as latest or stable.
const TAG: &str = r"(?::[A-Za-z0-9_][A-Za-z0-9._-]*|@sha256:[0-9a-f]{64})";
/// What may not stand right before an image, nor right after it.
const AHEAD: &str = r"(?:[^A-Za-z0-9_/-]|$)";

fn lead(c: char) -> bool {
    c.is_ascii_alphanumeric() || "_./:@$-".contains(c)
}

static URL: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[a-z][a-z0-9+.-]*://[^\s\x1c-\x1f]+").expect("a valid pattern"));
/// A path to a file, such as modules/local/vep/main.nf:23, is not an image.
static FILE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r"\.(nf|sh|md|py|config|ya?ml|json|txt|tsv|csv|vcf|gz|html|css|js|R|r):[0-9]+(-[0-9]+)?$",
    )
    .expect("a valid pattern")
});

fn bare_names(pins: &[Pin]) -> BTreeSet<String> {
    let mut bare: BTreeSet<String> = BASES.iter().map(|s| (*s).to_owned()).collect();
    for p in pins.iter().filter(|p| p.name.ends_with("_IMAGE")) {
        if let Some(i) = p.value.find([':', '@']) {
            let name = &p.value[..i];
            if !name.is_empty() && !name.contains('/') {
                bare.insert(name.to_owned());
            }
        }
    }
    bare
}

/// The two patterns: `ns/name:tag`, and an official name with a tag.
fn patterns(bare: &BTreeSet<String>) -> [Rule; 2] {
    let names: Vec<String> = bare.iter().map(|b| regex::escape(b)).collect();
    [
        Rule::new(
            &format!(r"([a-z][a-z0-9._-]*/[a-z][a-z0-9._/-]*{TAG}){AHEAD}"),
            Some(lead),
        )
        .expect("a valid pattern"),
        Rule::new(
            &format!(r"((?:{}){TAG}){AHEAD}", names.join("|")),
            Some(lead),
        )
        .expect("a valid pattern"),
    ]
}

/// The text of a file as Python's text mode reads it: strict UTF-8, any
/// newline turned into `\n`. None for a file that is not UTF-8 or not a
/// readable file.
fn read_text(path: &Path) -> Option<String> {
    let bytes = std::fs::read(path).ok()?;
    let text = String::from_utf8(bytes).ok()?;
    Some(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// The literal images of one file, with their line numbers; exempt
/// occurrences are taken off `left`.
fn scan_file(
    rel: &str,
    text: &str,
    rules: &[Rule; 2],
    left: &mut BTreeMap<(String, String), usize>,
    out: &mut dyn Write,
) -> bool {
    let mut bad = false;
    for (i, raw) in text.split('\n').enumerate() {
        // docker://image:tag (a step image in a workflow) is an image, not a URL.
        let line = raw.replace("docker://", " ");
        let line = URL.replace_all(&line, " ");
        for r in rules {
            for (s, e) in r.find_iter(&line) {
                // A full stop ends the sentence, not the tag.
                let image = line[s..e].trim_end_matches('.');
                if FILE.is_match(image) || image.starts_with("example/") {
                    continue; // example/ is the placeholder namespace of the self-tests
                }
                if let Some(n) = left.get_mut(&(rel.to_owned(), image.to_owned()))
                    && *n > 0
                {
                    *n -= 1;
                    continue;
                }
                let _ = writeln!(
                    out,
                    "FAIL: {rel}:{} names the image {image}; use its *_IMAGE variable from versions.env",
                    i + 1
                );
                bad = true;
            }
        }
    }
    bad
}

struct ScanArgs {
    root: String,
    exempt_files: Vec<String>,
    pairs: BTreeMap<(String, String), usize>,
}

fn parse_scan(args: &[String]) -> Result<ScanArgs, String> {
    let (mut root, mut files, mut pairs) = (None, None, None);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let v = it.next().ok_or_else(|| format!("{a} needs a value"));
        match a.as_str() {
            "--root" => root = Some(v?.clone()),
            "--exempt-files" => files = Some(v?.clone()),
            "--exempt-pairs" => pairs = Some(v?.clone()),
            s => return Err(format!("unknown argument {s}")),
        }
    }
    let exempt_files = files
        .ok_or("--exempt-files is required")?
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_owned))
        .collect();
    let mut counted = BTreeMap::new();
    for l in pairs.ok_or("--exempt-pairs is required")?.lines() {
        let w: Vec<&str> = l.split_whitespace().take(2).collect();
        if let [f, image] = w[..] {
            *counted.entry((f.to_owned(), image.to_owned())).or_insert(0) += 1;
        } else if !w.is_empty() {
            return Err(format!("an exempt pair needs a file and an image: {l}"));
        }
    }
    Ok(ScanArgs {
        root: root.ok_or("--root is required")?,
        exempt_files,
        pairs: counted,
    })
}

fn scan(a: &ScanArgs, list: &str, out: &mut dyn Write) -> Result<bool, String> {
    let files: Vec<&str> = list.lines().filter(|f| !f.is_empty()).collect();
    if files.is_empty() {
        let _ = writeln!(out, "FAIL: no file to scan");
        return Ok(false);
    }
    let root = Path::new(&a.root);
    let pins = versions_env::read(&root.join("versions.env"))?;
    let rules = patterns(&bare_names(&pins));
    let mut left = a.pairs.clone();
    let mut bad = false;
    for f in &files {
        let skipped = a
            .exempt_files
            .iter()
            .any(|s| *f == s || (s.ends_with('/') && f.starts_with(s.as_str())));
        if skipped {
            continue;
        }
        if let Some(text) = read_text(&root.join(f)) {
            bad |= scan_file(f, &text, &rules, &mut left, out);
        }
    }
    for ((f, image), n) in &left {
        if *n > 0 {
            let _ = writeln!(
                out,
                "FAIL: {f} no longer names {image}{}; drop its line from EXEMPT_PAIRS in scripts/ci/check-images.sh",
                if *n > 1 {
                    format!(" {n} more times")
                } else {
                    String::new()
                }
            );
            bad = true;
        }
    }
    let _ = writeln!(out, "scanned {} files", files.len());
    Ok(!bad)
}

/// `og dev scan images ARGS`.
pub fn cli_scan(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match parse_scan(args) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og dev scan images: error: {e}");
            return 2;
        }
    };
    let mut list = String::new();
    if let Err(e) = std::io::stdin().read_to_string(&mut list) {
        let _ = writeln!(err, "ERROR: cannot read the file list: {e}");
        return 2;
    }
    match scan(&a, &list, out) {
        Ok(true) => 0,
        Ok(false) => 1,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}

/// The release in an image's tag: `quay.io/x/samtools:1.24--h0_0` is 1.24.
fn release(value: &str) -> Option<&str> {
    let tag = value.rsplit('/').next().unwrap_or(value);
    let (_, t) = tag.split_once(':')?;
    Some(t.split("--").next().unwrap_or(t))
}

static HELD: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"\b(minimap2|samtools)=([0-9][0-9A-Za-z.]*)").expect("a valid pattern")
});

/// The check on parsed pins; prints one line per tool.
pub fn coupled(pins: &[Pin], out: &mut dyn Write) -> bool {
    let Some(image) = versions_env::get(pins, "MINIMAP2_SAMTOOLS_IMAGE") else {
        let _ = writeln!(
            out,
            "FAIL: MINIMAP2_SAMTOOLS_IMAGE is not set in versions.env"
        );
        return false;
    };
    let mut held: BTreeMap<&str, &str> = BTreeMap::new();
    for c in HELD.captures_iter(&image.rest) {
        let (tool, v) = (
            c.get(1).map_or("", |m| m.as_str()),
            c.get(2).map_or("", |m| m.as_str()),
        );
        held.insert(tool, v.trim_end_matches('.'));
    }
    let mut ok = true;
    for (tool, var) in [
        ("minimap2", "MINIMAP2_IMAGE"),
        ("samtools", "SAMTOOLS_IMAGE"),
    ] {
        let want = versions_env::get(pins, var).and_then(|p| release(&p.value));
        let got = held.get(tool).copied();
        match (want, got) {
            (Some(w), Some(g)) if !w.is_empty() && w == g => {
                let _ = writeln!(
                    out,
                    "OK: MINIMAP2_SAMTOOLS_IMAGE holds {tool} {g}, the release of {var}"
                );
            }
            _ => {
                let want = want.filter(|w| !w.is_empty());
                let _ = writeln!(
                    out,
                    "FAIL: MINIMAP2_SAMTOOLS_IMAGE says {tool}={} in its versions.env comment, but {var} is {}: pin a BioContainers mulled image with {tool} {} and name it in that comment (docs/updating.md, Held pins)",
                    got.unwrap_or("(nothing)"),
                    want.unwrap_or("(unset)"),
                    want.unwrap_or("?")
                );
                ok = false;
            }
        }
    }
    ok
}

/// `og dev versions coupled --env FILE`.
pub fn cli_coupled(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let [flag, path] = args else {
        let _ = write!(err, "{USAGE}");
        return 2;
    };
    if flag != "--env" {
        let _ = write!(err, "{USAGE}");
        return 2;
    }
    match versions_env::read(Path::new(path)) {
        Ok(pins) => u8::from(!coupled(&pins, out)),
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rules() -> [Rule; 2] {
        let pins = versions_env::parse("PYTHON_IMAGE=\"python:3.15.2\"\n").unwrap();
        patterns(&bare_names(&pins))
    }

    fn found(line: &str) -> Vec<String> {
        let mut left = BTreeMap::new();
        let mut out = Vec::new();
        scan_file("f", line, &rules(), &mut left, &mut out);
        String::from_utf8(out)
            .unwrap()
            .lines()
            .map(|l| {
                l.split(" names the image ")
                    .nth(1)
                    .unwrap()
                    .split(';')
                    .next()
                    .unwrap()
                    .to_owned()
            })
            .collect()
    }

    // Planted images are split by "~", which is dropped, so this file
    // names none of them whole.
    fn p(s: &str) -> String {
        s.replace('~', "")
    }

    #[test]
    fn literals_are_found_as_python_finds_them() {
        assert_eq!(
            found(&p("docker run staphb/bcftools:~1.21 x")),
            [p("staphb/bcftools:~1.21")]
        );
        assert_eq!(
            found(&p("runs staphb/samtools:~1.19.")),
            [p("staphb/samtools:~1.19")]
        );
        assert_eq!(found(&p("// container python:~3.12")), [p("python:~3.12")]);
        assert_eq!(found(&p("uses: docker://alpine:~3.8")), [p("alpine:~3.8")]);
        // A tag followed by a slash backs off to the part before a dot.
        assert_eq!(
            found(&p("x staphb/samtools:~1.2/x")),
            [p("staphb/samtools:~1")]
        );
        for clean in [
            "see https://quay.io/repository/biocontainers/samtools:1.0",
            "modules/local/vep/main.nf:23",
            "example/bcftools:1.0",
            "$VAR/name:tag",
            "python3 x.py",
        ] {
            assert!(found(clean).is_empty(), "{clean}: {:?}", found(clean));
        }
    }

    #[test]
    fn a_stale_exemption_fails() {
        let dir = std::env::temp_dir().join(format!("og-dev-images-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("versions.env"), "A_IMAGE=\"x/y:1\"\n").unwrap();
        std::fs::write(dir.join("a.md"), p("uses quay.io/x/y:~2\n")).unwrap();
        let a = ScanArgs {
            root: dir.to_string_lossy().into_owned(),
            exempt_files: vec![],
            pairs: BTreeMap::from([
                (("a.md".into(), p("quay.io/x/y:~2")), 1),
                (("a.md".into(), p("quay.io/x/y:~3")), 1),
            ]),
        };
        let mut out = Vec::new();
        let ok = scan(&a, "a.md\n", &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(!ok);
        assert!(
            out.contains(&p(
                "FAIL: a.md no longer names quay.io/x/y:~3; drop its line"
            )),
            "{out}"
        );
        assert!(!out.contains(&p("y:~2;")), "{out}");
    }

    #[test]
    fn coupled_releases() {
        let env = "MINIMAP2_IMAGE=\"q/minimap2:2.31--h1_0\"\nSAMTOOLS_IMAGE=\"q/samtools:1.24--h2_0\"\nMINIMAP2_SAMTOOLS_IMAGE=\"q/mulled:abc\"  # minimap2=2.31,samtools=1.24.\n";
        let mut out = Vec::new();
        assert!(coupled(&versions_env::parse(env).unwrap(), &mut out));
        let bumped = env.replace("samtools:1.24", "samtools:9.99");
        let mut out = Vec::new();
        assert!(!coupled(&versions_env::parse(&bumped).unwrap(), &mut out));
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("FAIL: MINIMAP2_SAMTOOLS_IMAGE says samtools=1.24 in its versions.env comment, but SAMTOOLS_IMAGE is 9.99:"), "{out}");
    }
}
