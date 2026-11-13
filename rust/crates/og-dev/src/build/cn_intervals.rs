//! `og dev build cn-intervals`: the interval rows of `bin/cn_loci.tsv`,
//! the copy-number checks step 35 runs with one `samtools bedcov` pass
//! (`og cn intervals`).
//!
//! The rows are:
//!
//! - every region row of ClinGen Dosage Sensitivity's recurrent CNV BED,
//!   version 2.1, GRCh38 (`CLINGEN_RECURRENT_CNV_URL` in `versions.env`,
//!   CC0). Its breakpoint and LCR rows (itemRgb `0,0,0`) are the segmental
//!   duplications that bound the regions: they are cut out of the targets
//!   and kept out of the controls. Its chrX rows are spelled `chrx` and are
//!   read as chrX; any other name that is not a GRCh38 chromosome stops the
//!   build;
//! - the rows the BED lacks (`HAND_ROWS`): two ClinGen curations (1p36,
//!   Xq28 MECP2), five carrier deletions at their ClinVar or HbVar GRCh38
//!   spans, and the Y chromosome's AZF regions with two single-copy Y
//!   controls.
//!
//! Targets are cut to pieces of at most 500 kb and lose every stretch inside
//! a breakpoint or LCR row, a GIAB v3.6 segmental duplication, an assembly
//! gap (GIAB's gaps with 15 kb on each side), a satellite or a region of a
//! model the pinned Parascopy data holds; the AZF targets keep their
//! segmental duplications, since the DAZ copies are one. Controls are two
//! 100 kb windows, one on each side of the cluster of ClinGen rows the
//! region belongs to, the nearest that stay clear of all of that and of
//! every target, within 3 Mb; a side with none (a chromosome end) is left
//! out.
//!
//! Every cut set is merged once per chromosome, so a subtraction or a
//! "clear of the cuts" test is a binary search.

use super::{Args, BuildError, CHROMS, Failure, chrom_index, each_line, input_error, int};
use og_io::Utf8;
use std::io::Write;
use std::path::Path;

/// The line after which the catalog's rows are this builder's.
pub const MARKER: &str =
    "# --- interval rows: written by scripts/ci/build-cn-intervals.sh; do not edit by hand ---";
/// How a marker line starts, whatever builder it names.
const MARKER_START: &[u8] = b"# --- interval rows: ";

const REGION_RGB: &str = "255,111,0";
const MASK_RGB: &str = "0,0,0";
/// GRCh38 chromosome lengths (the no-ALT analysis set's `.fai`), in
/// [`CHROMS`] order.
const LENGTHS: [i64; 24] = [
    248_956_422,
    242_193_529,
    198_295_559,
    190_214_555,
    181_538_259,
    170_805_979,
    159_345_973,
    145_138_636,
    138_394_717,
    133_797_422,
    135_086_622,
    133_275_309,
    114_364_328,
    107_043_718,
    101_991_189,
    90_338_345,
    83_257_441,
    80_373_285,
    58_617_616,
    64_444_167,
    46_709_983,
    50_818_468,
    156_040_895,
    57_227_415,
];
const PIECE: i64 = 500_000;
const FLANK: i64 = 100_000;
const GAP: i64 = 50_000;
const SEARCH: i64 = 3_000_000;
const STEP: i64 = 10_000;
const MIN_PIECE: i64 = 500;
const CHR_Y: usize = 23;

/// A 0-based, half-open stretch.
type Span = (i64, i64);

/// How a row's target depth is read.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Depth {
    /// Depth pooled over the targets.
    Pooled,
    /// The targets' depths added up (copies of one sequence).
    Sum,
}

impl Depth {
    fn name(self) -> &'static str {
        match self {
            Depth::Pooled => "pooled",
            Depth::Sum => "sum",
        }
    }
}

/// A row the BED lacks.
struct Hand {
    locus: &'static str,
    chrom: usize,
    /// (label, 0-based start, end); the first label names the row.
    targets: &'static [(&'static str, i64, i64)],
    /// (name, start, end, label), or `None` for flanks.
    controls: Option<&'static [(&'static str, i64, i64, &'static str)]>,
    depth: Depth,
    /// Whether segmental duplications are cut out of the targets.
    cut_segdups: bool,
    pmids: &'static str,
}

const Y_CONTROLS: &[(&str, i64, i64, &str)] = &[
    (
        "ZFY",
        2_935_200,
        2_982_579,
        "ZFY, single-copy control on Yp",
    ),
    (
        "NLGN4Y",
        14_522_572,
        14_845_654,
        "NLGN4Y, single-copy control on Yq",
    ),
];

/// The rows the ClinGen BED lacks. 1p36 (ISCA-37434) and Xq28 (ISCA-46304)
/// from ClinGen_region_curation_list_GRCh38.tsv as read on 2026-11-04;
/// del(GJB6-D13S1830), CFTRdele2,3 and the 621 bp HBB deletion at their
/// ClinVar GRCh38 spans; the 619 bp HBB deletion (HbVar 979, Orkin 1979)
/// and the Sicilian (delta-beta)0 deletion (HbVar 1035, Henthorn 1990) at
/// their HbVar spans, placed on NC_000011.10 by 50 bp on each side of the
/// NG_000007.3 breakpoints on 2026-11-04; the AZF regions by the single-copy
/// genes each deletion removes (AZFc as the four DAZ copies added up), gene
/// spans from Ensembl's REST lookup, checked on 2026-11-04. The alpha-globin
/// --SEA deletion is read by step 35's HBA table, not here.
const HAND_ROWS: &[Hand] = &[
    Hand {
        locus: "1p36_terminal",
        chrom: 0,
        targets: &[(
            "1p36 terminal region (includes GABRD), ClinGen ISCA-37434",
            898_702,
            6_229_913,
        )],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: ".",
    },
    Hand {
        locus: "Xq28_MECP2",
        chrom: 22,
        targets: &[(
            "Xq28 region (includes MECP2), ClinGen ISCA-46304: the MECP2 duplication",
            154_008_528,
            154_110_279,
        )],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: ".",
    },
    Hand {
        locus: "GJB6_D13S1830",
        chrom: 12,
        targets: &[("del(GJB6-D13S1830), ClinVar 5546", 20_223_037, 20_531_806)],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: "11807148",
    },
    Hand {
        locus: "CFTR_dele2_3",
        chrom: 6,
        targets: &[(
            "CFTRdele2,3 (21 kb), ClinVar 66105",
            117_498_312,
            117_519_392,
        )],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: ".",
    },
    Hand {
        locus: "HBB_exon3_del",
        chrom: 10,
        targets: &[(
            "HBB deletion of exon 3 (621 bp), ClinVar 3769040",
            5_225_254,
            5_225_875,
        )],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: ".",
    },
    Hand {
        locus: "HBB_619bp_del",
        chrom: 10,
        targets: &[(
            "HBB 619 bp deletion (South Asian), HbVar 979",
            5_225_388,
            5_226_007,
        )],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: "287080",
    },
    Hand {
        locus: "HBD_HBB_Sicilian_del",
        chrom: 10,
        targets: &[(
            "Sicilian (delta-beta)0 deletion (13.4 kb, HBD to HBB 3'), HbVar 1035",
            5_219_877,
            5_233_280,
        )],
        controls: None,
        depth: Depth::Pooled,
        cut_segdups: true,
        pmids: "2307466",
    },
    Hand {
        locus: "AZFa",
        chrom: CHR_Y,
        targets: &[
            ("AZFa (USP9Y and DDX3Y)", 12_537_649, 12_860_852),
            ("DDX3Y", 12_904_107, 12_920_478),
        ],
        controls: Some(Y_CONTROLS),
        depth: Depth::Pooled,
        cut_segdups: false,
        pmids: "37674303",
    },
    Hand {
        locus: "AZFb",
        chrom: CHR_Y,
        targets: &[
            ("AZFb (KDM5D and EIF1AY)", 19_702_297, 19_745_148),
            ("EIF1AY", 20_575_693, 20_593_154),
        ],
        controls: Some(Y_CONTROLS),
        depth: Depth::Pooled,
        cut_segdups: false,
        pmids: "37674303",
    },
    Hand {
        locus: "AZFc_DAZ",
        chrom: CHR_Y,
        targets: &[
            (
                "AZFc (DAZ1 to DAZ4 added up; gr/gr and b2/b3 leave two)",
                23_129_354,
                23_199_010,
            ),
            ("DAZ2", 23_219_446, 23_291_356),
            ("DAZ3", 24_763_068, 24_813_492),
            ("DAZ4", 24_833_842, 24_907_040),
        ],
        controls: Some(Y_CONTROLS),
        depth: Depth::Sum,
        cut_segdups: false,
        pmids: "37674303",
    },
];

/// A region row of the ClinGen BED.
#[derive(Debug)]
struct Region {
    chrom: usize,
    start: i64,
    end: i64,
    name: String,
}

/// The stretches cut out of targets and controls, per chromosome.
#[derive(Default, Clone)]
struct Cuts {
    /// ClinGen breakpoint and LCR rows.
    mask: Vec<Span>,
    segdup: Vec<Span>,
    gap: Vec<Span>,
    satellite: Vec<Span>,
    parascopy: Vec<Span>,
}

/// The downloads, read.
pub struct Inputs {
    regions: Vec<Region>,
    cuts: Vec<Cuts>,
}

/// `chrx` and `chry` (the ClinGen BED's spelling) as chrX and chrY; any
/// other name that is not a GRCh38 chromosome is an error.
fn norm_chrom(c: &str) -> Result<usize, BuildError> {
    let fixed = match c {
        "chrx" => "chrX",
        "chry" => "chrY",
        _ => c,
    };
    chrom_index(fixed).ok_or_else(|| {
        BuildError::Input(format!(
            "chromosome {} is not a GRCh38 chromosome name",
            py_repr(c)
        ))
    })
}

fn py_repr(s: &str) -> String {
    pycompat::value::py_repr(&pycompat::json::Value::Str(s.to_owned()))
        .unwrap_or_else(|_| format!("'{s}'"))
}

/// The region rows and breakpoint rows of the ClinGen BED's lines.
fn read_clingen<'a>(
    lines: impl IntoIterator<Item = &'a str>,
    cuts: &mut [Cuts],
) -> Result<Vec<Region>, BuildError> {
    let mut regions = Vec::new();
    for line in lines {
        let strip = pycompat::text::strip;
        if strip(line, None).is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 9 {
            return Err(BuildError::Input(format!(
                "ClinGen BED: a row with {} columns: {}",
                c.len(),
                py_repr(strip(line, None))
            )));
        }
        let chrom = norm_chrom(c[0])?;
        let (start, end) = (
            int(c[1], "ClinGen BED start")?,
            int(c[2], "ClinGen BED end")?,
        );
        let (name, rgb) = (strip(c[3], None), strip(c[8], None));
        if end <= start {
            return Err(BuildError::Input(format!(
                "ClinGen BED: {name}: empty span"
            )));
        }
        match rgb {
            REGION_RGB => regions.push(Region {
                chrom,
                start,
                end,
                name: name.to_owned(),
            }),
            MASK_RGB => cuts[chrom].mask.push((start, end)),
            _ => {
                return Err(BuildError::Input(format!(
                    "ClinGen BED: {name}: itemRgb {} is neither a region nor a breakpoint row",
                    py_repr(rgb)
                )));
            }
        }
    }
    if regions.is_empty() {
        return Err(BuildError::Input("ClinGen BED: no region row".to_owned()));
    }
    Ok(regions)
}

/// The rows of a BED on chr1-22, X and Y, per chromosome.
fn read_bed(
    path: &Path,
    put: impl Fn(&mut Cuts) -> &mut Vec<Span>,
    cuts: &mut [Cuts],
) -> Result<(), BuildError> {
    each_line(path, Utf8::Strict, |line| {
        if line.starts_with('#') || pycompat::text::strip(line, None).is_empty() {
            return Ok(());
        }
        let mut c = line.split('\t');
        let (seq, s, e) = (c.next(), c.next(), c.next());
        let (Some(seq), Some(s), Some(e)) = (seq, s, e) else {
            return Err(BuildError::Input(format!(
                "{}: a row with under 3 columns",
                path.display()
            )));
        };
        if let Some(i) = chrom_index(seq) {
            let span = (int(s, "BED start")?, int(e, "BED end")?);
            put(&mut cuts[i]).push(span);
        }
        Ok(())
    })
}

/// The regions of Parascopy's models: `NAME CHROM:START-END ...`, 1-based.
fn read_parascopy(path: &Path, cuts: &mut [Cuts]) -> Result<(), BuildError> {
    each_line(path, Utf8::Strict, |line| {
        if line.starts_with('#') || pycompat::text::strip(line, None).is_empty() {
            return Ok(());
        }
        let bad = || {
            BuildError::Input(format!(
                "{}: not NAME CHROM:START-END DETAIL: {line:?}",
                path.display()
            ))
        };
        let f = pycompat::text::split(line, None, None);
        let [_, reg, _] = f[..] else {
            return Err(bad());
        };
        let (c, span) = reg.split_once(':').ok_or_else(bad)?;
        let (s, e) = span.split_once('-').ok_or_else(bad)?;
        if let Some(i) = chrom_index(c) {
            cuts[i]
                .parascopy
                .push((int(s, "Parascopy start")? - 1, int(e, "Parascopy end")?));
        }
        Ok(())
    })
}

/// Sorted, with overlapping and touching spans joined.
fn merge(mut spans: Vec<Span>) -> Vec<Span> {
    spans.sort_unstable();
    let mut out: Vec<Span> = Vec::with_capacity(spans.len());
    for (s, e) in spans {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// The parts of `[s, e)` outside the merged `cuts`.
fn subtract(s: i64, e: i64, cuts: &[Span]) -> Vec<Span> {
    let mut out = Vec::new();
    let mut pos = s;
    let first = cuts.partition_point(|c| c.1 <= s);
    for &(cs, ce) in &cuts[first..] {
        if cs >= e {
            break;
        }
        if cs > pos {
            out.push((pos, cs));
        }
        pos = pos.max(ce);
    }
    if pos < e {
        out.push((pos, e));
    }
    out
}

/// Each span of `MIN_PIECE` bp or more cut into equal pieces of at most
/// `PIECE` bp.
fn pieces(spans: &[Span]) -> Vec<Span> {
    let mut out = Vec::new();
    for &(s, e) in spans {
        let len = e - s;
        if len < MIN_PIECE {
            continue;
        }
        // Ceiling divisions of positive numbers.
        let n = (len + PIECE - 1) / PIECE;
        let step = (len + n - 1) / n;
        // PLANTED, reverted in the next commit: each piece ends one base short.
        out.extend((0..n).map(|i| (s + i * step, e.min(s + (i + 1) * step) - 1)));
    }
    out
}

/// Whether `[s, e)` touches none of the merged `cuts`.
fn clear(s: i64, e: i64, cuts: &[Span]) -> bool {
    let i = cuts.partition_point(|c| c.1 <= s);
    cuts.get(i).is_none_or(|c| c.0 >= e)
}

/// The nearest `FLANK` window `GAP` beyond `edge` on the left (`left`) or
/// right, clear of `cuts`, within `SEARCH`.
fn flank(chrom: usize, edge: i64, left: bool, cuts: &[Span]) -> Option<Span> {
    for k in 0..=SEARCH / STEP {
        let (s, e) = if left {
            let e = edge - GAP - k * STEP;
            (e - FLANK, e)
        } else {
            let s = edge + GAP + k * STEP;
            (s, s + FLANK)
        };
        if s < 0 || e > LENGTHS[chrom] {
            return None;
        }
        if clear(s, e, cuts) {
            return Some((s, e));
        }
    }
    None
}

/// The row id of a ClinGen region name: its band, then the words of its
/// parentheses that are not an "includes" note.
fn locus_id(name: &str) -> String {
    let band = pycompat::text::split(name, None, None)
        .first()
        .copied()
        .unwrap_or_default();
    let b = name.as_bytes();
    let mut groups: Vec<&str> = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'(' {
            let j = b[i + 1..]
                .iter()
                .position(|&c| c == b'(' || c == b')')
                .map(|k| i + 1 + k);
            if let Some(j) = j
                && b[j] == b')'
            {
                let g = &name[i + 1..j];
                if !pycompat::text::strip(g, None).starts_with("includes") {
                    groups.push(g);
                }
                i = j + 1;
                continue;
            }
        }
        i += 1;
    }
    let mut words = String::new();
    let mut in_run = false;
    for ch in groups.join(" ").chars() {
        if ch.is_ascii_alphanumeric() || ch == '.' || ch == '-' {
            words.push(ch);
            in_run = false;
        } else if !in_run {
            words.push('_');
            in_run = true;
        }
    }
    let words = words.trim_matches('_');
    if words.is_empty() {
        band.to_owned()
    } else {
        format!("{band}_{words}")
    }
}

/// A name made safe for a target or control cell.
fn label(text: &str) -> String {
    let t: String = text
        .chars()
        .map(|c| {
            if matches!(c, ';' | '\t' | '\n') {
                ','
            } else {
                c
            }
        })
        .collect();
    pycompat::text::strip(&t, None).to_owned()
}

/// A control window: (name, start, end, label).
type Control = (String, i64, i64, String);

/// The catalog rows, as Python wrote them, in their order.
struct Builder<'a> {
    inputs: &'a Inputs,
    /// Per chromosome, the ClinGen rows that touch or overlap, merged
    /// (each row's end taken one base further).
    clusters: Vec<Vec<Span>>,
    /// Per chromosome, every cut (with and without segmental duplications)
    /// and the cuts a control keeps out of.
    cuts_all: Vec<Vec<Span>>,
    cuts_no_segdups: Vec<Vec<Span>>,
    flank_cuts: Vec<Vec<Span>>,
    rows: Vec<String>,
    seen: std::collections::HashMap<String, usize>,
}

impl<'a> Builder<'a> {
    fn new(inputs: &'a Inputs) -> Builder<'a> {
        let n = CHROMS.len();
        let mut clusters = vec![Vec::new(); n];
        let mut flank_cuts = vec![Vec::new(); n];
        for r in &inputs.regions {
            clusters[r.chrom].push((r.start, r.end + 1));
            flank_cuts[r.chrom].push((r.start, r.end));
        }
        for h in HAND_ROWS {
            flank_cuts[h.chrom].extend(h.targets.iter().map(|&(_, s, e)| (s, e)));
        }
        let mut cuts_all = Vec::with_capacity(n);
        let mut cuts_no_segdups = Vec::with_capacity(n);
        for (i, c) in inputs.cuts.iter().enumerate() {
            clusters[i].extend(c.mask.iter().map(|&(s, e)| (s, e + 1)));
            let all: Vec<Span> = [&c.mask, &c.segdup, &c.gap, &c.satellite, &c.parascopy]
                .into_iter()
                .flatten()
                .copied()
                .collect();
            flank_cuts[i].extend(&all);
            cuts_all.push(merge(all));
            cuts_no_segdups.push(merge(
                [&c.mask, &c.gap, &c.parascopy]
                    .into_iter()
                    .flatten()
                    .copied()
                    .collect(),
            ));
        }
        Builder {
            inputs,
            clusters: clusters.into_iter().map(merge).collect(),
            cuts_all,
            cuts_no_segdups,
            flank_cuts: flank_cuts.into_iter().map(merge).collect(),
            rows: Vec::new(),
            seen: std::collections::HashMap::new(),
        }
    }

    /// The two flank controls of the cluster that holds `[s, e]`.
    fn flanks(&self, chrom: usize, s: i64, e: i64) -> Result<Vec<Control>, BuildError> {
        let (cs, ce) = self.clusters[chrom]
            .iter()
            .copied()
            .find(|&(a, b)| a <= s && e <= b)
            .unwrap_or((s, e));
        let mut out = Vec::new();
        for (left, edge, n) in [(true, cs, "L"), (false, ce, "R")] {
            if let Some((a, b)) = flank(chrom, edge, left, &self.flank_cuts[chrom]) {
                let side = if left { "left" } else { "right" };
                out.push((n.to_owned(), a, b, format!("flank {side} of the region")));
            }
        }
        if out.is_empty() {
            return Err(BuildError::Input(format!(
                "{}:{s}-{e}: no flank clear of the masks within {} Mb",
                CHROMS[chrom],
                SEARCH / 1_000_000
            )));
        }
        Ok(out)
    }

    fn add(
        &mut self,
        locus: &str,
        chrom: usize,
        targets: &[(i64, i64, String)],
        controls: &[Control],
        depth: Depth,
        pmids: &str,
    ) -> Result<(), BuildError> {
        if targets.is_empty() {
            return Err(BuildError::Input(format!(
                "{locus}: no target stretch is left outside the masks"
            )));
        }
        let n = self.seen.entry(locus.to_owned()).or_insert(0);
        *n += 1;
        let locus = if *n > 1 {
            format!("{locus}_{n}")
        } else {
            locus.to_owned()
        };
        let tg: Vec<String> = targets
            .iter()
            .enumerate()
            .map(|(i, (s, e, lab))| {
                let tail = if lab.is_empty() {
                    String::new()
                } else {
                    format!(":{lab}")
                };
                format!("{locus}_t{}:{s}-{e}{tail}", i + 1)
            })
            .collect();
        let ct: Vec<String> = controls
            .iter()
            .map(|(n, s, e, lab)| format!("{locus}_{n}:{s}-{e}:{lab}"))
            .collect();
        let y = chrom == CHR_Y;
        let cells = [
            locus.as_str(),
            "interval",
            CHROMS[chrom],
            &tg.join(";"),
            &ct.join(";"),
            ".",
            depth.name(),
            if y { "pooled" } else { "flanks" },
            "2",
            "0",
            ".",
            ".",
            ".",
            if y { "4" } else { "8" },
            if y { "." } else { "1.5" },
            "round:0.3",
            ".",
            ".",
            "bedcov",
            "cn-depth",
            "cnvpytor:cnvs",
            pmids,
        ];
        self.rows.push(cells.join("\t"));
        Ok(())
    }

    /// Targets with the first piece labelled.
    fn labelled(segs: &[Span], name: &str) -> Vec<(i64, i64, String)> {
        segs.iter()
            .enumerate()
            .map(|(i, &(a, b))| (a, b, if i == 0 { label(name) } else { String::new() }))
            .collect()
    }

    fn build(mut self) -> Result<Vec<String>, BuildError> {
        let mut regions: Vec<&Region> = self.inputs.regions.iter().collect();
        regions.sort_by_key(|r| (r.chrom, r.start, r.end));
        for r in regions {
            let segs = pieces(&subtract(r.start, r.end, &self.cuts_all[r.chrom]));
            let targets = Self::labelled(&segs, &r.name);
            let controls = self.flanks(r.chrom, r.start, r.end)?;
            self.add(
                &locus_id(&r.name),
                r.chrom,
                &targets,
                &controls,
                Depth::Pooled,
                ".",
            )?;
        }
        for h in HAND_ROWS {
            let cuts = if h.cut_segdups {
                &self.cuts_all[h.chrom]
            } else {
                &self.cuts_no_segdups[h.chrom]
            };
            let mut targets = Vec::new();
            for &(name, s, e) in h.targets {
                let segs = match h.depth {
                    Depth::Pooled => pieces(&subtract(s, e, cuts)),
                    Depth::Sum => vec![(s, e)],
                };
                if h.depth == Depth::Sum
                    && !clear(s, e, &merge(self.inputs.cuts[h.chrom].parascopy.clone()))
                {
                    return Err(BuildError::Input(format!(
                        "{}: {name} overlaps a Parascopy model",
                        h.locus
                    )));
                }
                targets.extend(Self::labelled(&segs, name));
            }
            let controls = match h.controls {
                None => self.flanks(h.chrom, h.targets[0].1, h.targets[h.targets.len() - 1].2)?,
                Some(c) => c
                    .iter()
                    .enumerate()
                    .map(|(i, &(_, s, e, lab))| (format!("c{}", i + 1), s, e, lab.to_owned()))
                    .collect(),
            };
            self.add(h.locus, h.chrom, &targets, &controls, h.depth, h.pmids)?;
        }
        Ok(self.rows)
    }
}

impl Inputs {
    /// The downloads: the ClinGen BED's lines, the three GIAB BEDs and the
    /// Parascopy model list.
    pub fn read(
        clingen: &Path,
        segdups: &Path,
        gaps: &Path,
        satellites: &Path,
        parascopy: &Path,
    ) -> Result<Inputs, BuildError> {
        let mut cuts = vec![Cuts::default(); CHROMS.len()];
        let lines: Vec<String> = super::lines(clingen, Utf8::Replace)?
            .collect::<Result<_, _>>()
            .map_err(|e| input_error(clingen, &e))?;
        let regions = read_clingen(lines.iter().map(String::as_str), &mut cuts)?;
        read_bed(segdups, |c| &mut c.segdup, &mut cuts)?;
        read_bed(gaps, |c| &mut c.gap, &mut cuts)?;
        read_bed(satellites, |c| &mut c.satellite, &mut cuts)?;
        read_parascopy(parascopy, &mut cuts)?;
        Ok(Inputs { regions, cuts })
    }

    /// The interval rows, without line endings, marker excluded.
    pub fn rows(&self) -> Result<Vec<String>, BuildError> {
        Builder::new(self).build()
    }
}

/// The catalog with its rows from the marker line on replaced: the lines
/// before a marker line are kept (all of them when there is none).
fn splice(catalog: &[u8], rows: &[String]) -> Vec<u8> {
    let mut out = Vec::with_capacity(catalog.len());
    for line in catalog.split_inclusive(|&c| c == b'\n') {
        if line.starts_with(MARKER_START) && line.ends_with(b"\n") {
            break;
        }
        out.extend_from_slice(line);
    }
    for line in std::iter::once(MARKER).chain(rows.iter().map(String::as_str)) {
        out.extend_from_slice(line.as_bytes());
        out.push(b'\n');
    }
    out
}

/// How two catalogs differ: (lines only in `old`, lines only in `new`),
/// each line counted as often as it appears.
fn line_diff(old: &[u8], new: &[u8]) -> (usize, usize) {
    let mut count: std::collections::HashMap<&[u8], i64> = std::collections::HashMap::new();
    for l in old.split_inclusive(|&c| c == b'\n') {
        *count.entry(l).or_default() += 1;
    }
    for l in new.split_inclusive(|&c| c == b'\n') {
        *count.entry(l).or_default() -= 1;
    }
    let only_old = count
        .values()
        .filter(|&&v| v > 0)
        .map(|&v| v as usize)
        .sum();
    let only_new = count
        .values()
        .filter(|&&v| v < 0)
        .map(|&v| (-v) as usize)
        .sum();
    (only_old, only_new)
}

pub(super) fn cli(
    args: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<u8, Failure> {
    let a = Args::parse(
        args,
        &[
            "--clingen",
            "--segdups",
            "--gaps",
            "--satellites",
            "--parascopy",
            "--catalog",
        ],
        &["--check", "--stdout"],
    )?;
    let catalog = a.get("--catalog");
    if catalog.is_some() == a.has("--stdout") {
        return Err(Failure::Usage(
            "give one of --catalog FILE and --stdout".to_owned(),
        ));
    }
    if a.has("--check") && catalog.is_none() {
        return Err(Failure::Usage("--check needs --catalog".to_owned()));
    }
    let p = |f: &str| a.need(f).map(Path::new);
    let inputs = Inputs::read(
        p("--clingen")?,
        p("--segdups")?,
        p("--gaps")?,
        p("--satellites")?,
        p("--parascopy")?,
    )?;
    let rows = inputs.rows()?;
    let Some(catalog) = catalog else {
        let mut text = String::new();
        for line in std::iter::once(MARKER).chain(rows.iter().map(String::as_str)) {
            text.push_str(line);
            text.push('\n');
        }
        let _ = out.write_all(text.as_bytes());
        return Ok(0);
    };
    let path = Path::new(catalog);
    let old = std::fs::read(path).map_err(|e| input_error(path, &e))?;
    let new = splice(&old, &rows);
    if a.has("--check") {
        if new != old {
            let (only_old, only_new) = line_diff(&old, &new);
            let how = if only_old == 0 && only_new == 0 {
                "the same lines in another order".to_owned()
            } else {
                format!("{only_old} line(s) only in the file, {only_new} only in the rebuild")
            };
            let _ = writeln!(
                err,
                "{catalog} is not what the pinned downloads give ({how}); rerun scripts/ci/build-cn-intervals.sh"
            );
            return Ok(1);
        }
        let _ = writeln!(
            out,
            "{catalog}: the {} interval rows are what the pinned downloads give",
            rows.len()
        );
        return Ok(0);
    }
    super::write_atomic(path, &new)?;
    let _ = writeln!(out, "{catalog}: {} interval rows written", rows.len());
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small ClinGen BED: the 22q11.2 A-D region between its breakpoint
    /// rows, and the Xp22.31 region with its breakpoint, spelled `chrx`.
    const SMALL_BED: &str = "track name=\"test\"\n#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\nchr22\t18159724\t18924717\tBP-A\t0\t.\t18159724\t18924717\t0,0,0\nchr22\t18924718\t21111383\t22q11.2 recurrent (DGS) region (proximal, A-D) (includes TBX1)\t0\t.\t18924718\t21111383\t255,111,0\nchr22\t21111384\t21562827\tBP-D\t0\t.\t21111384\t21562827\t0,0,0\nchrx\t6527561\t6537770\tXp22.31 distal breakpoint\t0\t.\t6527561\t6537770\t0,0,0\nchrx\t6537771\t8156913\tXp22.31 recurrent region (includes STS)\t0\t.\t6537771\t8156913\t255,111,0\n";

    fn small(bed: &str) -> Result<Vec<String>, BuildError> {
        let mut cuts = vec![Cuts::default(); CHROMS.len()];
        let regions = read_clingen(bed.lines(), &mut cuts)?;
        Inputs { regions, cuts }.rows()
    }

    fn cell(row: &str, i: usize) -> &str {
        row.split('\t').nth(i).unwrap()
    }

    /// (start, end) of each target of a row.
    fn spans(row: &str) -> Vec<Span> {
        cell(row, 3)
            .split(';')
            .map(|t| {
                let (s, e) = t.split(':').nth(1).unwrap().split_once('-').unwrap();
                (s.parse().unwrap(), e.parse().unwrap())
            })
            .collect()
    }

    #[test]
    fn the_small_bed_gives_its_two_regions_then_the_hand_rows() {
        let rows = small(SMALL_BED).unwrap();
        let ids: Vec<(&str, &str)> = rows.iter().map(|r| (cell(r, 0), cell(r, 2))).collect();
        assert_eq!(
            ids[..2],
            [("22q11.2_DGS_proximal_A-D", "chr22"), ("Xp22.31", "chrX")]
        );
        let hand: Vec<&str> = ids[2..].iter().map(|x| x.0).collect();
        let want: Vec<&str> = HAND_ROWS.iter().map(|h| h.locus).collect();
        assert_eq!(hand, want, "the hand rows follow, in their order");
        assert!(rows.iter().all(|r| r.split('\t').count() == 22));
        assert!(rows.iter().all(|r| cell(r, 20) == "cnvpytor:cnvs"));
    }

    #[test]
    fn an_unknown_contig_stops_the_build() {
        let e = small(&SMALL_BED.replace("chrx", "chrQ")).unwrap_err();
        assert!(
            e.to_string().contains("'chrQ' is not a GRCh38 chromosome"),
            "{e}"
        );
        assert_eq!(e.code(), 2, "exit 2, as the Python builder did");
    }

    #[test]
    fn no_target_piece_overlaps_a_breakpoint_row() {
        let rows = small(SMALL_BED).unwrap();
        for (s, e) in spans(&rows[0]) {
            assert!(!(s < 18_924_717 && 18_159_724 < e), "{s}-{e} on BP-A");
            assert!(!(s < 21_562_827 && 21_111_384 < e), "{s}-{e} on BP-D");
        }
        // 2.19 Mb in pieces of at most 500 kb, the first labelled.
        let t = cell(&rows[0], 3);
        assert_eq!(t.split(';').count(), 5);
        assert!(t.starts_with("22q11.2_DGS_proximal_A-D_t1:18924718-"));
        assert!(
            t.split(';')
                .next()
                .unwrap()
                .ends_with(":22q11.2 recurrent (DGS) region (proximal, A-D) (includes TBX1)")
        );
    }

    #[test]
    fn controls_flank_the_cluster_and_never_sit_on_a_target() {
        let rows = small(SMALL_BED).unwrap();
        // 22q11.2: the cluster runs from BP-A to BP-D, so the controls lie
        // beyond those, 50 kb out.
        assert_eq!(
            cell(&rows[0], 4),
            "22q11.2_DGS_proximal_A-D_L:18009724-18109724:flank left of the region;22q11.2_DGS_proximal_A-D_R:21612828-21712828:flank right of the region"
        );
        for r in &rows {
            if cell(r, 2) == "chrY" {
                continue;
            }
            for c in cell(r, 4).split(';') {
                let (s, e) = c.split(':').nth(1).unwrap().split_once('-').unwrap();
                let (cs, ce): (i64, i64) = (s.parse().unwrap(), e.parse().unwrap());
                for o in &rows {
                    if cell(o, 2) == cell(r, 2) {
                        assert!(
                            spans(o).iter().all(|&(ts, te)| !(ts < ce && cs < te)),
                            "{} control on {}",
                            cell(r, 0),
                            cell(o, 0)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn hand_rows_keep_their_spans_and_columns() {
        let rows = small(SMALL_BED).unwrap();
        let get = |l: &str| rows.iter().find(|r| cell(r, 0) == l).unwrap().clone();
        assert_eq!(spans(&get("HBB_619bp_del")), [(5_225_388, 5_226_007)]);
        assert_eq!(
            spans(&get("HBD_HBB_Sicilian_del")),
            [(5_219_877, 5_233_280)]
        );
        let az = get("AZFc_DAZ");
        assert_eq!(
            (cell(&az, 6), cell(&az, 7), cell(&az, 13), cell(&az, 14)),
            ("sum", "pooled", "4", ".")
        );
        assert_eq!(spans(&az).len(), 4, "the DAZ copies are never cut");
        assert!(cell(&az, 4).starts_with("AZFc_DAZ_c1:2935200-2982579:ZFY"));
        let x = get("Xq28_MECP2");
        assert_eq!(
            (cell(&x, 7), cell(&x, 13), cell(&x, 14)),
            ("flanks", "8", "1.5")
        );
    }

    #[test]
    fn a_parascopy_model_is_cut_out_of_a_target() {
        let mut cuts = vec![Cuts::default(); CHROMS.len()];
        let regions = read_clingen(SMALL_BED.lines(), &mut cuts).unwrap();
        cuts[21].parascopy.push((19_000_000, 19_100_000));
        let rows = Inputs { regions, cuts }.rows().unwrap();
        assert!(
            spans(&rows[0])
                .iter()
                .all(|&(s, e)| e <= 19_000_000 || s >= 19_100_000)
        );
        assert!(spans(&rows[0]).iter().any(|&(_, e)| e == 19_000_000));
    }

    #[test]
    fn a_parascopy_model_on_a_daz_copy_stops_the_build() {
        let mut cuts = vec![Cuts::default(); CHROMS.len()];
        let regions = read_clingen(SMALL_BED.lines(), &mut cuts).unwrap();
        cuts[CHR_Y].parascopy.push((24_800_000, 24_800_100));
        let e = Inputs { regions, cuts }.rows().unwrap_err();
        assert!(
            e.to_string()
                .contains("AZFc_DAZ: DAZ3 overlaps a Parascopy model"),
            "{e}"
        );
    }

    #[test]
    fn bad_rows_are_refused() {
        let mut c = vec![Cuts::default(); CHROMS.len()];
        assert!(
            read_clingen(["chr1\t1\t2"], &mut c)
                .unwrap_err()
                .to_string()
                .contains("a row with 3 columns")
        );
        let e = read_clingen(["chr1\t5\t5\tx\t0\t.\t5\t5\t255,111,0"], &mut c).unwrap_err();
        assert!(e.to_string().contains("x: empty span"));
        let e = read_clingen(["chr1\t1\t5\tx\t0\t.\t1\t5\t1,2,3"], &mut c).unwrap_err();
        assert!(e.to_string().contains("itemRgb '1,2,3'"));
        let e = read_clingen(["chr1\t1\t5\tx\t0\t.\t1\t5\t0,0,0"], &mut c).unwrap_err();
        assert!(e.to_string().contains("no region row"));
    }

    #[test]
    fn interval_helpers() {
        assert_eq!(
            merge(vec![(5, 8), (1, 3), (3, 4), (7, 10)]),
            [(1, 4), (5, 10)]
        );
        assert_eq!(
            subtract(0, 20, &[(2, 4), (6, 8), (30, 40)]),
            [(0, 2), (4, 6), (8, 20)]
        );
        assert_eq!(subtract(3, 7, &[(0, 4), (6, 9)]), [(4, 6)]);
        assert!(subtract(3, 7, &[(0, 9)]).is_empty());
        assert_eq!(
            pieces(&[(0, 499), (0, 1_000_001)]),
            [(0, 333_334), (333_334, 666_668), (666_668, 1_000_001)]
        );
        assert_eq!(pieces(&[(0, 500)]), [(0, 500)]);
        assert!(clear(4, 6, &[(0, 4), (6, 9)]));
        assert!(!clear(3, 6, &[(0, 4), (6, 9)]));
        assert_eq!(flank(0, 1_000_000, true, &[]), Some((850_000, 950_000)));
        assert_eq!(
            flank(0, 1_000_000, true, &[(900_000, 945_000)]),
            Some((800_000, 900_000))
        );
        assert_eq!(
            flank(0, 100_000, true, &[]),
            None,
            "past the chromosome start"
        );
    }

    #[test]
    fn locus_ids_and_labels() {
        assert_eq!(
            locus_id("22q11.2 recurrent (DGS) region (proximal, A-D) (includes TBX1)"),
            "22q11.2_DGS_proximal_A-D"
        );
        assert_eq!(
            locus_id("Xp22.31 recurrent region (includes STS)"),
            "Xp22.31"
        );
        assert_eq!(locus_id("1q21.1 (a (b) c)"), "1q21.1_b");
        assert_eq!(label(" a;b\tc "), "a,b,c");
    }

    #[test]
    fn a_repeated_locus_gets_a_number() {
        let bed = "chr1\t10000000\t10100000\t1p (X) one\t0\t.\t0\t0\t255,111,0\nchr1\t20000000\t20100000\t1p (X) two\t0\t.\t0\t0\t255,111,0\n";
        let rows = small(bed).unwrap();
        assert_eq!((cell(&rows[0], 0), cell(&rows[1], 0)), ("1p_X", "1p_X_2"));
        assert!(cell(&rows[1], 3).starts_with("1p_X_2_t1:"));
    }

    #[test]
    fn splice_keeps_the_head_and_check_counts_both_sides() {
        let rows = vec!["a\t1".to_owned(), "b\t2".to_owned()];
        let old = b"# head\nx\n# --- interval rows: written by an older builder ---\nold\n";
        let new = splice(old, &rows);
        assert_eq!(
            new,
            format!("# head\nx\n{MARKER}\na\t1\nb\t2\n").into_bytes()
        );
        assert_eq!(
            splice(b"x\n", &rows),
            format!("x\n{MARKER}\na\t1\nb\t2\n").into_bytes()
        );
        let added = [new.as_slice(), b"c\t3\n"].concat();
        assert_eq!(line_diff(&new, &added), (0, 1), "an added line is counted");
        assert_eq!(line_diff(&added, &new), (1, 0));
        let swapped = format!("# head\nx\n{MARKER}\nb\t2\na\t1\n").into_bytes();
        assert_eq!(line_diff(&new, &swapped), (0, 0));
    }
}
