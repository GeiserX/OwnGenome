//! Each step's `metric<TAB>value` table and stdout line for a block, and
//! the line of the interval table. The HBA and LPA_KIV2 layouts are those
//! `bin/hba_copy_number.py` and `bin/lpa_kiv2.py` wrote; a `call psv` row
//! gets the PSV layout and every other row the generic one.
//!
//! Each layout matches on the block's status and reason exhaustively. The
//! HBA layout wrote `status called` with no count for a reason it did not
//! name (a catalog edit could give `high_mapq0_share` or `below_floor`), and
//! the LPA layout failed on one; both now say `not assessed: <reason>`.

use crate::catalog::Locus;
use crate::engine::{Block, Reason, Status};
use crate::{d, f2, fixed};

/// The rows of a step table and its one stdout line.
pub type Table = (Vec<(String, String)>, String);

/// How a row's table is laid out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    Generic,
    Hba,
    LpaKiv2,
    Psv,
}

impl Layout {
    /// The layout of `row`.
    pub fn of(row: &Locus) -> Layout {
        if !row.psv.is_empty() {
            Layout::Psv
        } else {
            match row.id.as_str() {
                "HBA" => Layout::Hba,
                "LPA_KIV2" => Layout::LpaKiv2,
                _ => Layout::Generic,
            }
        }
    }
}

/// The step table of `row` for `b`.
pub fn table(row: &Locus, b: &Block) -> Table {
    match Layout::of(row) {
        Layout::Psv => psv(row, b),
        Layout::Hba => hba(row, b),
        Layout::LpaKiv2 => lpa(row, b),
        Layout::Generic => generic(row, b),
    }
}

fn kv(rows: Vec<(&str, String)>) -> Vec<(String, String)> {
    rows.into_iter().map(|(k, v)| (k.to_owned(), v)).collect()
}

fn reason(b: &Block) -> String {
    b.reason.map_or_else(|| ".".to_owned(), |r| r.to_string())
}

fn or_dot(s: &str) -> String {
    if s.is_empty() { ".".into() } else { s.into() }
}

fn generic(row: &Locus, b: &Block) -> Table {
    let rows = kv(vec![
        ("status", b.status.to_string()),
        ("copies", d(b.copies)),
        ("estimate", f2(b.estimate)),
        ("target_depth", f2(b.target_depth)),
        ("control_depth", f2(b.control_depth)),
        ("reason", reason(b)),
        ("reading", or_dot(&b.reading)),
    ]);
    let why = b
        .reason
        .map_or_else(|| b.reading.clone(), |r| r.to_string());
    (rows, format!("{}: {} ({why})", row.id, b.status))
}

/// `"<reason>"` for a reason a layout has no words for.
fn bare(b: &Block) -> String {
    format!("not assessed: {}", reason(b))
}

fn hba(row: &Locus, b: &Block) -> Table {
    const ORDER: [&str; 9] = [
        "status",
        "alpha_globin_genes",
        "estimate",
        "hba2_depth",
        "hba1_depth",
        "flank_left_depth",
        "flank_right_depth",
        "mapq0_share",
        "message",
    ];
    if b.reason == Some(Reason::NoData) {
        let msg = "not assessed: the BAM has no chr16 reads, so HBA1 and HBA2 have none";
        let rows = ORDER
            .iter()
            .map(|k| {
                let v = match *k {
                    "status" => "not assessed",
                    "message" => msg,
                    _ => ".",
                };
                ((*k).to_owned(), v.to_owned())
            })
            .collect();
        return (rows, format!("HBA1/HBA2: {msg}"));
    }
    let depth = |n: &str| b.depths.get(n).copied();
    let r2 = |x: Option<f64>| {
        x.map_or_else(
            || ".".to_owned(),
            |v| pycompat::fmt::repr(pycompat::fmt::round(v, 2).unwrap_or(v)),
        )
    };
    let (left, right) = (depth("flank_left"), depth("flank_right"));
    let one = |x: Option<f64>| x.map_or_else(|| ".".to_owned(), |v| fixed(v, 1));
    let share = b.mapq0_share;
    let est = b.estimate;
    let (status, genes, message) = match (b.status, b.reason) {
        (Status::Called, _) => ("called", d(b.copies), b.reading.clone()),
        (_, Some(Reason::LowControlDepth)) => (
            "not assessed",
            ".".into(),
            format!(
                "not assessed: flank depth under {}x ({}x and {}x), too low to compare",
                pycompat::fmt::general(row.depth().min_control, 6),
                one(left),
                one(right)
            ),
        ),
        (_, Some(Reason::ControlsDisagree)) => (
            "not assessed",
            ".".into(),
            format!(
                "not assessed: the two flanks disagree ({}x and {}x); a deletion larger than --SEA, or uneven coverage, would skew the count",
                one(left),
                one(right)
            ),
        ),
        (_, Some(Reason::LowMapq0Share)) => (
            "not assessed",
            ".".into(),
            format!(
                "not assessed: reads of mapping quality 0 make up {}% of the depth over HBA1 and HBA2, so the aligner left out the reads that fit both genes (as novoalign does by default) and the count would come out low; realign with minimap2 (step 02) or BWA-MEM",
                one(share.map(|s| 100.0 * s))
            ),
        ),
        (_, Some(Reason::BetweenCounts)) => {
            let e = est.expect("an estimate between counts");
            let n = e.trunc() as i64;
            (
                "not assessed",
                ".".into(),
                format!(
                    "not assessed: the estimate {} lies between {n} and {} alpha-globin genes",
                    fixed(e, 2),
                    n + 1
                ),
            )
        }
        _ => ("not assessed", ".".into(), bare(b)),
    };
    let values = [
        status.to_owned(),
        genes,
        f2(est),
        r2(depth("HBA2")),
        r2(depth("HBA1")),
        r2(left),
        r2(right),
        share.map_or_else(|| ".".to_owned(), |s| fixed(s, 3)),
        message,
    ];
    let line = format!(
        "HBA1/HBA2: {} (estimate {}; HBA2 {}x, HBA1 {}x, flanks {}x and {}x; MAPQ 0 share over the genes {})",
        values[8], values[2], values[3], values[4], values[5], values[6], values[7]
    );
    let rows = ORDER
        .iter()
        .zip(values)
        .map(|(k, v)| ((*k).to_owned(), v))
        .collect();
    (rows, line)
}

/// The advice line of step 46's table.
pub const LPA_ADVICE: &str = "A research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 European Atherosclerosis Society consensus advises measuring it at least once in every adult.";

fn lpa(row: &Locus, b: &Block) -> Table {
    let mut status = "not assessed".to_owned();
    let (mut copies, mut lo, mut hi, mut class) = (
        ".".to_owned(),
        ".".to_owned(),
        ".".to_owned(),
        ".".to_owned(),
    );
    let (mut kd, mut fd, mut sh) = (".".to_owned(), ".".to_owned(), ".".to_owned());
    let message;
    if b.reason == Some(Reason::NoData) {
        message = "not assessed: the BAM has no chr6, so LPA has no reads".to_owned();
    } else {
        let k = b.target_depth;
        let f = b.control_depth;
        let one = |x: Option<f64>| x.map_or_else(|| ".".to_owned(), |v| fixed(v, 1));
        kd = f2(k);
        fd = f2(f);
        sh = b
            .mapq0_share
            .map_or_else(|| ".".to_owned(), |s| fixed(s, 3));
        message = match (b.status, b.reason) {
            (Status::Called, _) => match (b.estimate, b.ci95) {
                (Some(est), Some((l, h))) => {
                    status = "estimated".into();
                    copies = fixed(est, 1);
                    lo = fixed(l, 1);
                    hi = fixed(h, 1);
                    class = b.class.map_or(".", |c| c.as_str()).to_owned();
                    format!(
                        "about {} KIV-2 repeats in total (95% interval {lo} to {hi}): {}",
                        fixed(est, 0),
                        b.reading
                    )
                }
                _ => bare(b),
            },
            (_, Some(Reason::LowControlDepth)) => {
                status = "too little depth".into();
                format!(
                    "too little depth: {}x over the single-copy LPA exons, under the {}x the count needs",
                    one(f),
                    pycompat::fmt::general(row.depth().min_control, 6)
                )
            }
            (_, Some(Reason::LowMapq0Share)) => format!(
                "not assessed: reads of mapping quality 0 make up {}% of the depth over the KIV-2 exons, so the aligner left out the reads that fit several KIV-2 copies (as novoalign does by default) and the count would come out low; realign with minimap2 (step 02) or BWA-MEM",
                one(b.mapq0_share.map(|s| 100.0 * s))
            ),
            (_, Some(Reason::BelowFloor)) => format!(
                "not assessed: the KIV-2 exons have {}x against {}x over the single-copy exons, too little for any repeat count (a deletion or a mapping problem)",
                one(k),
                one(f)
            ),
            _ => bare(b),
        };
    }
    let line = format!(
        "LPA KIV-2: {message} (KIV-2 exons {kd}x, single-copy exons {fd}x, MAPQ 0 share {sh})"
    );
    let rows = kv(vec![
        ("status", status),
        ("kiv2_copies", copies),
        ("ci95_low", lo),
        ("ci95_high", hi),
        ("class", class),
        ("kiv2_depth", kd),
        ("flank_depth", fd),
        ("mapq0_share", sh),
        ("message", message),
        ("advice", LPA_ADVICE.to_owned()),
    ]);
    (rows, line)
}

/// `.` for a called PSV row, else `not assessed: ` and why; a low MAPQ 0
/// share is spelled out, because the BAM, not the genome, stopped the call.
fn psv_message(row: &Locus, b: &Block) -> String {
    if b.status == Status::Called {
        return ".".into();
    }
    if b.reason == Some(Reason::LowMapq0Share) {
        let targets: Vec<&str> = row
            .targets
            .iter()
            .map(|s| s.name.strip_suffix("_BODY").unwrap_or(&s.name))
            .collect();
        let share = b.mapq0_share.unwrap_or(f64::NAN);
        let min = row.depth().mapq0_min.unwrap_or(f64::NAN);
        return format!(
            "not assessed: the BAM has few or no multi-mapping reads here (reads of mapping quality 0 make up {}% of the depth over {}, under {}%), a possible sign of an aligner that leaves out the reads that fit both copies (as novoalign does by default); the module total may then come out low and read as a deletion; realign with minimap2 (step 02) or BWA-MEM",
            fixed(100.0 * share, 1),
            targets.join(" and "),
            pycompat::fmt::general(100.0 * min, 6)
        );
    }
    format!("not assessed: {}", reason(b))
}

fn psv(row: &Locus, b: &Block) -> Table {
    let mut rows = kv(vec![
        ("status", b.status.to_string()),
        ("cyp21a2_copies", d(b.copies)),
        ("module_total", d(b.module_total)),
        ("expected_total", d(b.expected)),
        ("estimate", f2(b.estimate)),
        ("finding", b.finding.map_or(".", |f| f.as_str()).to_owned()),
        ("conversion_sites", or_dot(&b.conversion_sites.join(","))),
    ]);
    for s in &b.sites {
        rows.push((
            format!("psv_{}", s.name),
            format!(
                "{} gene, {} pseudogene, {} other reads; copies {} ({}); ratio {}",
                s.gene_reads,
                s.pseudo_reads,
                s.other_reads,
                d(s.copies),
                f2(s.estimate),
                f2(s.ratio)
            ),
        ));
    }
    rows.extend(kv(vec![
        ("target_depth", f2(b.target_depth)),
        ("control_depth", f2(b.control_depth)),
        ("reason", reason(b)),
        ("reading", or_dot(&b.reading)),
        (
            "mapq0_share",
            b.mapq0_share
                .map_or_else(|| ".".to_owned(), |s| fixed(s, 3)),
        ),
        ("message", psv_message(row, b)),
    ]));
    let why = match (b.reason, b.finding) {
        (Some(r), _) => r.to_string(),
        (None, Some(f)) => f.as_str().to_owned(),
        (None, None) => "None".to_owned(),
    };
    let line = format!(
        "{}: {} ({why}; {} gene copies, module total {} of {})",
        row.id,
        b.status,
        d(b.copies),
        d(b.module_total),
        d(b.expected)
    );
    (rows, line)
}

/// The interval table's columns.
pub const INTERVAL_COLUMNS: [&str; 11] = [
    "locus", "name", "chrom", "start", "end", "status", "copies", "expected", "change", "estimate",
    "reason",
];

/// One line of the interval table.
pub fn interval_line(row: &Locus, b: &Block) -> Vec<String> {
    let tg = &row.targets;
    let change = match (b.copies, b.expected) {
        (Some(c), Some(e)) if c < e => "deletion",
        (Some(c), Some(e)) if c > e => "duplication",
        (Some(_), Some(_)) => "no change",
        _ => ".",
    };
    let name = if tg[0].label.is_empty() {
        row.id.clone()
    } else {
        tg[0].label.clone()
    };
    vec![
        row.id.clone(),
        name,
        row.chrom.clone(),
        tg.iter()
            .map(|s| s.start)
            .min()
            .expect("a target")
            .to_string(),
        tg.iter()
            .map(|s| s.end)
            .max()
            .expect("a target")
            .to_string(),
        b.status.to_string(),
        d(b.copies),
        d(b.expected),
        change.to_owned(),
        f2(b.estimate),
        reason(b),
    ]
}
