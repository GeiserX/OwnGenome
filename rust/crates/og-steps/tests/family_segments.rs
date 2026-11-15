//! `og family segments` (step 38's chromosome map) and the map in `og
//! family report`, on a synthetic three-generation family: two
//! grandparents, their son (the father), the mother and a grandchild.
//! Genotypes are drawn from known haplotypes, so the truth is known.
//! It replaces `tests/test_family_segments.py`.
//!
//! 1. trio (father and mother sequenced, the father's parents too): on
//!    chr1 the grandchild got the grandfather's copy up to a crossover
//!    planted at 20 Mb and the grandmother's after it. The map has those two
//!    paternal segments, each boundary within 50 kb of the crossover; three
//!    planted genotype errors inside the first segment stay inside it
//!    (discordant 3), not new segments; chr2 has no crossover, and its 10 Mb
//!    gap with no sites splits it into two segments from the grandmother;
//!    the mother's side, with no grandparents sequenced, is one stretch per
//!    block; a son's chrX comes from his mother only, outside the PARs;
//! 2. duo (the mother named, not sequenced): the paternal crossover is
//!    still found, and the maternal side is the copy that did not come from
//!    the father, over the same stretches;
//! 3. the DNA Painter CSV has DNA Painter's own header, a row per segment,
//!    chromosome names without `chr`, and the grandparent as Ancestor;
//! 4. the SVG is well-formed XML with a bar per segment;
//! 5. without `--par-bed` chrX is left out and the autosomes are mapped;
//! 6. `og family report` writes the files and a "Chromosome map" section in
//!    family_report.md that counts the crossover; a second run writes the
//!    same bytes.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

const CROSSOVER: i64 = 20_000_000;
const STEP: i64 = 2_000;
const ERRORS: [i64; 3] = [5_000_000, 9_000_000, 13_000_000];

fn par_bed() -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../assets/par_grch38.bed")
        .display()
        .to_string()
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-family-seg-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn family(args: &[&str]) -> u8 {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    og_steps::family::cli(&args)
}

/// A small deterministic generator (xorshift64*), so the family is the
/// same on every run.
struct Rng(u64);

impl Rng {
    fn bit(&mut self) -> u8 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        u8::from(self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 63 == 1)
    }
}

fn gt(a: u8, b: u8) -> String {
    format!("{}/{}", a.min(b), a.max(b))
}

fn members_tsv(path: &Path, mum_sequenced: bool) -> String {
    let rows = [
        ("kid", "male", "dad", "mum", ""),
        ("dad", "male", "gpa", "gma", "mum"),
        ("mum", "female", "", "", "dad"),
        ("gpa", "male", "", "", "gma"),
        ("gma", "female", "", "", "gpa"),
    ];
    let mut t = String::from("sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n");
    for (s, sex, fa, mo, pa) in rows {
        let seq = if s == "mum" && !mum_sequenced {
            "no"
        } else {
            "yes"
        };
        t += &format!("{s}\t{sex}\t{fa}\t{mo}\t{pa}\tyes\t{seq}\n");
    }
    fs::write(path, t).unwrap();
    path.display().to_string()
}

/// Sites on chr1 (a crossover in the father's meiosis at CROSSOVER, and
/// the three planted errors), chr2 (no crossover, a 10 Mb gap with no
/// sites) and chrX (a son's single X from his mother).
fn family_vcf(path: &Path) -> String {
    let samples = ["gpa", "gma", "dad", "mum", "kid"];
    let mut rng = Rng(38);
    let blocks: [(&str, &[(i64, i64)]); 3] = [
        ("chr1", &[(1_000_000, 41_000_000)]),
        ("chr2", &[(1_000_000, 11_000_000), (21_000_000, 31_000_000)]),
        ("chrX", &[(5_000_000, 15_000_000)]),
    ];
    let mut out = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n\
         ##contig=<ID=chr2,length=242193529>\n##contig=<ID=chrX,length=156040895>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for s in samples {
        out += &format!("\t{s}");
    }
    out.push('\n');
    for (chrom, spans) in blocks {
        for &(lo, hi) in spans {
            let mut pos = lo;
            while pos < hi {
                let mut h = [0u8; 6];
                for x in &mut h {
                    *x = rng.bit();
                }
                let [mut gpa1, mut gpa2, mut gma1, mut gma2, mut mum1, mut mum2] = h;
                let planted = chrom == "chr1" && ERRORS.contains(&pos);
                if planted {
                    // a site that tells everything: the father 0/1 from gpa
                    // 0/0 and gma 1/1, the mother 0/0
                    (gpa1, gpa2, gma1, gma2, mum1, mum2) = (0, 0, 1, 1, 0, 0);
                }
                // the father carries gpa1 (from the grandfather) and gma1
                let g: [String; 5] = if chrom == "chrX" {
                    [
                        gpa1.to_string(),
                        gt(gma1, gma2),
                        gma1.to_string(),
                        gt(mum1, mum2),
                        mum1.to_string(),
                    ]
                } else {
                    let from_dad = if planted {
                        gma1 // a genotype error that looks like the other copy
                    } else if chrom == "chr1" && pos < CROSSOVER {
                        gpa1
                    } else {
                        gma1
                    };
                    [
                        gt(gpa1, gpa2),
                        gt(gma1, gma2),
                        gt(gpa1, gma1),
                        gt(mum1, mum2),
                        gt(from_dad, mum1),
                    ]
                };
                out += &format!(
                    "{chrom}\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{}\n",
                    g.join("\t")
                );
                pos += STEP;
            }
        }
    }
    fs::write(path, out).unwrap();
    path.display().to_string()
}

type Row = std::collections::HashMap<String, String>;

fn read_tsv(path: &Path) -> Vec<Row> {
    let text = fs::read_to_string(path).unwrap_or_default();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

fn n(r: &Row, k: &str) -> i64 {
    r[k].parse().unwrap()
}

fn near(x: i64, target: i64) -> bool {
    (x - target).abs() <= 50_000
}

fn pick<'a>(rows: &'a [Row], side: &str, chrom: &str) -> Vec<&'a Row> {
    rows.iter()
        .filter(|r| r["side"] == side && r["chrom"] == chrom)
        .collect()
}

/// The tags of `xml` balance (each open tag has its close, in order) and
/// it has one root element; returns the `<rect class="seg"` bars.
fn svg_bars(xml: &str) -> Result<usize, String> {
    let mut stack: Vec<&str> = Vec::new();
    let mut roots = 0;
    let mut bars = 0;
    let mut rest = xml;
    while let Some(i) = rest.find('<') {
        let j = rest[i..].find('>').ok_or("a tag is not closed")? + i;
        let tag = &rest[i + 1..j];
        rest = &rest[j + 1..];
        if tag.starts_with('?') || tag.starts_with('!') {
            continue;
        }
        if let Some(name) = tag.strip_prefix('/') {
            let open = stack.pop().ok_or(format!("</{name}> closes nothing"))?;
            if open != name.trim() {
                return Err(format!("<{open}> closed by </{name}>"));
            }
            continue;
        }
        let name = tag.split([' ', '\n', '\t', '/']).next().unwrap_or_default();
        if stack.is_empty() {
            roots += 1;
        }
        if name == "rect" && tag.contains("class=\"seg\"") {
            bars += 1;
        }
        if !tag.ends_with('/') {
            stack.push(name);
        }
    }
    if !stack.is_empty() || roots != 1 {
        return Err(format!("unclosed {stack:?}, {roots} root elements"));
    }
    Ok(bars)
}

#[test]
fn the_svg_check_can_fail() {
    assert_eq!(svg_bars("<svg><rect class=\"seg\"/></svg>"), Ok(1));
    assert!(svg_bars("<svg><g></svg>").is_err());
    assert!(svg_bars("<svg/><svg/>").is_err());
}

#[test]
fn trio_duo_no_par_and_report() {
    let d = scratch("all");
    let vcf = family_vcf(&d.join("family.vcf"));
    let members = members_tsv(&d.join("members.tsv"), true);
    let par = par_bed();

    // 1. trio
    let out = d.join("trio");
    let rc = family(&[
        "segments",
        "--members",
        &members,
        "--vcf",
        &vcf,
        "--par-bed",
        &par,
        "--out-dir",
        out.to_str().unwrap(),
    ]);
    assert_eq!(rc, 0, "og family segments exits 0");
    let rows = read_tsv(&out.join("kid_segments.tsv"));
    assert!(!rows.is_empty(), "kid_segments.tsv written");
    let pat1 = pick(&rows, "paternal", "chr1");
    assert_eq!(pat1.len(), 2, "chr1, paternal: two segments: {pat1:?}");
    let (a, b) = (pat1[0], pat1[1]);
    assert_eq!(
        (a["ancestor"].as_str(), b["ancestor"].as_str()),
        ("gpa", "gma"),
        "chr1: first from the grandfather, then from the grandmother"
    );
    assert!(
        near(n(a, "end"), CROSSOVER) && n(a, "end") < CROSSOVER,
        "chr1: the first segment ends within 50 kb before the crossover: {a:?}"
    );
    assert!(
        near(n(b, "start"), CROSSOVER) && n(b, "start") >= CROSSOVER,
        "chr1: the second starts within 50 kb after it: {b:?}"
    );
    assert!(
        near(n(a, "start"), 1_000_000) && near(n(b, "end"), 41_000_000),
        "chr1: the segments span the chromosome's sites"
    );
    assert_eq!(
        (a["discordant"].as_str(), b["discordant"].as_str()),
        ("2", "0"),
        "chr1: the three planted errors are discordant sites of the first segment"
    );
    assert_eq!(a["basis"], "parent and grandparents");
    let pat2 = pick(&rows, "paternal", "chr2");
    assert_eq!(
        pat2.len(),
        2,
        "chr2, paternal: exactly two segments: {pat2:?}"
    );
    for (r, (s, e)) in pat2
        .iter()
        .zip([(1_000_000, 11_000_000), (21_000_000, 31_000_000)])
    {
        assert!(
            r["ancestor"] == "gma" && near(n(r, "start"), s) && near(n(r, "end"), e),
            "chr2: split by the gap, both from the grandmother: {r:?}"
        );
    }
    assert!(
        pick(&rows, "paternal", "chrX").is_empty(),
        "chrX: nothing from the father of a son"
    );
    let matx = pick(&rows, "maternal", "chrX");
    assert!(
        matx.len() == 1 && matx[0]["parent"] == "mum",
        "chrX: one maternal segment, from the mother: {matx:?}"
    );
    let mat: Vec<&Row> = rows
        .iter()
        .filter(|r| r["side"] == "maternal" && (r["chrom"] == "chr1" || r["chrom"] == "chr2"))
        .collect();
    assert!(
        mat.len() == 3
            && mat
                .iter()
                .all(|r| r["basis"] == "parent" && r["ancestor"] == "."),
        "maternal side, no grandparents sequenced: one stretch per block: {mat:?}"
    );
    assert!(
        out.join("dad_segments.tsv").exists(),
        "dad has a map of his own"
    );
    assert!(
        !out.join("mum_segments.tsv").exists(),
        "no map for a member without a sequenced parent"
    );

    // 3. the DNA Painter CSV
    let csv = fs::read_to_string(out.join("kid_dnapainter.csv")).unwrap();
    let lines = csv.split_inclusive('\n').map(Ok::<_, std::io::Error>);
    let r: Vec<Vec<String>> =
        pycompat::csv::Reader::new(lines, ',', pycompat::csv::Quoting::Minimal)
            .map(Result::unwrap)
            .collect();
    assert_eq!(
        r[0],
        [
            "Chr",
            "Start",
            "End",
            "cM",
            "SNPs",
            "Match",
            "Confidence",
            "Group",
            "Side",
            "Notes",
            "Colour",
            "Ancestor"
        ],
        "DNA Painter CSV: DNA Painter's own header"
    );
    let body = &r[1..];
    assert_eq!(body.len(), rows.len(), "DNA Painter CSV: a row per segment");
    assert!(
        body.iter().all(|x| ["1", "2", "X"].contains(&x[0].as_str())
            && ["maternal", "paternal"].contains(&x[8].as_str())),
        "DNA Painter CSV: chromosomes without 'chr', sides maternal or paternal"
    );
    let first1 = body.iter().find(|x| x[0] == "1").unwrap();
    assert!(
        first1[7] == "gpa" && first1[11] == "1",
        "the first chr1 segment is in the grandparent's group, a known ancestor: {first1:?}"
    );
    assert!(
        body.iter().any(|x| x[7] == "gma" && x[11] == "1"),
        "the other grandparent has a group of its own"
    );
    assert!(
        body.iter().all(|x| x[11] == "1" || x[11].is_empty()),
        "Ancestor is only 1 or empty"
    );

    // 4. the SVG
    let svg = fs::read_to_string(out.join("kid_segments.svg")).unwrap();
    assert_eq!(
        svg_bars(&svg),
        Ok(rows.len()),
        "SVG: well-formed, a bar per segment"
    );

    // 2. duo
    let duo_members = members_tsv(&d.join("members_duo.tsv"), false);
    let out = d.join("duo");
    let rc = family(&[
        "segments",
        "--members",
        &duo_members,
        "--vcf",
        &vcf,
        "--par-bed",
        &par,
        "--out-dir",
        out.to_str().unwrap(),
    ]);
    assert_eq!(rc, 0, "duo: exits 0");
    let rows = read_tsv(&out.join("kid_segments.tsv"));
    let pat1 = pick(&rows, "paternal", "chr1");
    let anc: Vec<&str> = pat1.iter().map(|r| r["ancestor"].as_str()).collect();
    assert_eq!(anc, ["gpa", "gma"], "duo, chr1: the crossover is found");
    assert!(
        near(n(pat1[0], "end"), CROSSOVER) && near(n(pat1[1], "start"), CROSSOVER),
        "duo, chr1: boundary within 50 kb of the crossover: {pat1:?}"
    );
    let mat: Vec<&Row> = rows.iter().filter(|r| r["side"] == "maternal").collect();
    let chroms: BTreeSet<&str> = mat.iter().map(|r| r["chrom"].as_str()).collect();
    assert!(
        !mat.is_empty()
            && mat
                .iter()
                .all(|r| r["basis"] == "other copy" && r["parent"] == "mum")
            && chroms == BTreeSet::from(["chr1", "chr2"]),
        "duo: the maternal side is the other copy, on chr1 and chr2 only: {mat:?}"
    );
    let want = [
        ("chr1", 1_000_000, 41_000_000),
        ("chr2", 1_000_000, 11_000_000),
        ("chr2", 21_000_000, 31_000_000),
    ];
    assert!(
        mat.len() == 3
            && mat.iter().zip(want).all(|(r, (c, s, e))| r["chrom"] == c
                && near(n(r, "start"), s)
                && near(n(r, "end"), e)),
        "duo: the other copy covers chr1 and the two chr2 stretches: {mat:?}"
    );

    // 5. no --par-bed
    let out = d.join("nopar");
    let rc = family(&[
        "segments",
        "--members",
        &members,
        "--vcf",
        &vcf,
        "--out-dir",
        out.to_str().unwrap(),
    ]);
    assert_eq!(rc, 0, "no --par-bed: exits 0");
    let rows = read_tsv(&out.join("kid_segments.tsv"));
    assert!(
        !rows.is_empty() && rows.iter().all(|r| r["chrom"] != "chrX"),
        "no --par-bed: chrX is left out, the autosomes are still mapped"
    );

    // 6. og family report, twice
    let rel = d.join("relationships.tsv");
    fs::write(
        &rel,
        "sample_a\tsample_b\tdeclared\texpected_relatedness\trelatedness\tibs0_rate\tsites\tverdict\treason\n",
    )
    .unwrap();
    let outs = [d.join("report1"), d.join("report2")];
    for out in &outs {
        let rc = family(&[
            "report",
            "--members",
            &members,
            "--relationships",
            rel.to_str().unwrap(),
            "--vcf",
            &vcf,
            "--par-bed",
            &par,
            "--out-dir",
            out.to_str().unwrap(),
        ]);
        assert_eq!(rc, 0, "og family report exits 0");
    }
    let md = fs::read_to_string(outs[0].join("family_report.md")).unwrap();
    assert!(md.contains("## Chromosome map"), "a Chromosome map section");
    let line = md
        .lines()
        .find(|x| x.starts_with("- kid: from the father (dad)"))
        .unwrap_or_default();
    // The first "; N crossover" or "; N crossovers" followed by a comma.
    let count = line.match_indices("; ").find_map(|(i, _)| {
        let rest = &line[i + 2..];
        let digits = rest.len() - rest.trim_start_matches(|c: char| c.is_ascii_digit()).len();
        let tail = &rest[digits..];
        (digits > 0 && (tail.starts_with(" crossover,") || tail.starts_with(" crossovers,")))
            .then(|| &rest[..digits])
    });
    assert_eq!(
        count,
        Some("1"),
        "the section counts kid's one crossover: {line}"
    );
    for name in [
        "kid_segments.tsv",
        "kid_dnapainter.csv",
        "kid_segments.svg",
        "dad_segments.tsv",
    ] {
        let a = fs::read(outs[0].join(name)).unwrap();
        let b = fs::read(outs[1].join(name)).unwrap();
        assert!(
            a == b,
            "report writes {name}, the same bytes on a second run"
        );
    }
    fs::remove_dir_all(d).unwrap();
}
