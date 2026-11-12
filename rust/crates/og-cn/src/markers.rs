//! Copy-aware markers: two well-sourced markers sit in sequence the genome
//! holds more than once, so a gVCF genotype at one position cannot read
//! them. Each is read from the BAM by the step that counts its copies: the
//! reads at every GRCh38 copy of the site are added up, and the share that
//! carries the marker's base is a fraction of copies, not a diploid
//! genotype.
//!
//! - `LPA_KIV2_4925`: the LPA KIV-2 4925G>A splice variant (Coassin et al.
//!   2017, PMID 28444229), read by step 46 beside the KIV-2 count.
//! - `SMN1_27134`: SMN1 g.27134T>G (rs143838139; Luo et al. 2014, PMID
//!   23788250), read by step 35 beside the SMN1 psCN.
//!
//! The values are those `bin/traits_report.py` held (MARKERS), each copy
//! checked against the 21-base context and the dbSNP record by
//! `marker_errors` (a unit test here, and `tests/test_traits_report.py`
//! compares the records with what collect_summary reads).

use crate::input::each_pileup_site;
use crate::{Error, fixed};

/// One copy-aware marker.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Marker {
    pub key: &'static str,
    pub name: &'static str,
    pub chrom: &'static str,
    /// The plus-strand base read at each copy, and the marker's base.
    pub reference: u8,
    pub alt: u8,
    /// dbSNP b157's record of one copy: rsID, chrom, POS, REF, ALT.
    pub dbsnp: (&'static str, &'static str, i64, &'static str, &'static str),
    /// Every GRCh38 copy of the site, 1-based.
    pub copies: &'static [i64],
    /// The GRCh38 sequence of 10 bases on each side of every copy.
    pub context: &'static str,
    /// Fewer reads over all copies together: too little depth.
    pub min_depth: u64,
    /// Reads with the marker's base for 'present' ...
    pub min_alt: u64,
    /// ... and a share of the reads at least this.
    pub min_fraction: f64,
    /// This many reads with the marker's base or fewer: 'not seen'.
    pub max_absent: u64,
}

/// The LPA KIV-2 4925G>A splice variant: the last base of the second KIV-2
/// exon, G>A on LPA's minus strand, so C>T on the plus strand; one copy at
/// the end of each of the six KIV-2 units of GRCh38 and one at the
/// near-identical last exon of KIV-1. dbSNP names the first KIV-2 copy.
pub const LPA_KIV2_4925: Marker = Marker {
    key: "LPA_KIV2_4925",
    name: "LPA KIV-2 4925G>A",
    chrom: "chr6",
    reference: b'C',
    alt: b'T',
    dbsnp: ("rs1484919544", "chr6", 160_640_667, "C", "T"),
    copies: &[
        160_646_214,
        160_640_667,
        160_635_123,
        160_629_574,
        160_624_030,
        160_618_484,
        160_612_932,
    ],
    context: "GACTCCTTACCTTGTTCGGAA",
    min_depth: 100,
    min_alt: 5,
    min_fraction: 0.008,
    max_absent: 2,
};

/// SMN1 g.27134T>G (NM_000344.4:c.*3+80T>G, intron 7, plus strand); SMN2
/// holds the same 21 bases 875,420 bases upstream.
pub const SMN1_27134: Marker = Marker {
    key: "SMN1_27134",
    name: "SMN1 g.27134T>G",
    chrom: "chr5",
    reference: b'T',
    alt: b'G',
    dbsnp: ("rs143838139", "chr5", 70_952_074, "T", "G"),
    copies: &[70_952_074, 70_076_654],
    context: "TTTTTGAACATTTAAAAAGTT",
    min_depth: 10,
    min_alt: 3,
    min_fraction: 0.05,
    max_absent: 1,
};

/// Every marker, by key.
pub const MARKERS: [Marker; 2] = [LPA_KIV2_4925, SMN1_27134];

/// The marker named `key`.
pub fn marker(key: &str) -> Option<&'static Marker> {
    MARKERS.iter().find(|m| m.key == key)
}

/// A marker's reading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Present,
    NotSeen,
    Uncertain,
    TooLittleDepth,
    NotAssessed,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Present => "present",
            Status::NotSeen => "not seen",
            Status::Uncertain => "uncertain",
            Status::TooLittleDepth => "too little depth",
            Status::NotAssessed => "not assessed",
        }
    }
}

/// Reads with the reference base, with the marker's base, and others.
pub type Counts = (u64, u64, u64);

/// A marker's reading and the numbers behind it, as the table prints them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reading {
    pub status: Status,
    pub alt_reads: String,
    pub depth: String,
    /// `alt / depth` at 4 decimals, or `.`.
    pub fraction: String,
}

/// The copy sites as samtools regions, `chrom:pos-pos`, in the marker's order.
pub fn regions(m: &Marker) -> Vec<String> {
    m.copies
        .iter()
        .map(|p| format!("{}:{p}-{p}", m.chrom))
        .collect()
}

/// The reads over every copy of the marker from `samtools mpileup` lines
/// at its sites, without a reference: a read that carries the marker's base
/// counts wherever the aligner put it.
pub fn read_pileup(path: &str, m: &Marker) -> Result<Counts, Error> {
    let mut n = [0u64; 4];
    each_pileup_site(path, m.chrom, &|p| m.copies.contains(&p), |_, c| {
        for (a, b) in n.iter_mut().zip(c) {
            *a += b;
        }
    })?;
    let at = |b: u8| n["ACGT".find(b as char).expect("a base")];
    let (r, a) = (at(m.reference), at(m.alt));
    Ok((r, a, n.iter().sum::<u64>() - r - a))
}

/// The marker's reading from its counts (`None`: the pileup could not be
/// read). The depth counts the reads with the reference or the marker's base.
pub fn assess(m: &Marker, counts: Option<Counts>) -> Reading {
    let Some((r, a, _)) = counts else {
        return Reading {
            status: Status::NotAssessed,
            alt_reads: ".".into(),
            depth: ".".into(),
            fraction: ".".into(),
        };
    };
    let depth = r + a;
    let frac = (depth > 0).then(|| a as f64 / depth as f64);
    let status = if depth < m.min_depth {
        Status::TooLittleDepth
    } else if a >= m.min_alt && frac.is_some_and(|f| f >= m.min_fraction) {
        Status::Present
    } else if a <= m.max_absent {
        Status::NotSeen
    } else {
        Status::Uncertain
    };
    Reading {
        status,
        alt_reads: a.to_string(),
        depth: depth.to_string(),
        fraction: frac.map_or_else(|| ".".into(), |f| fixed(f, 4)),
    }
}

/// The printed fraction times `copies`, at one decimal (the reading of a
/// present or uncertain marker; `.` otherwise).
fn copies_of(a: &Reading, copies: Option<f64>) -> String {
    match (a.status, copies, pycompat::num::parse_float(&a.fraction)) {
        (Status::Present | Status::Uncertain, Some(c), Ok(f)) => fixed(f * c, 1),
        _ => ".".into(),
    }
}

const NULL_LIMIT: &str = "A short-read fraction cannot say which chromosome or which repeat carries it, so it cannot say whether it sits on the small or the large isoform.";

/// The KIV-2 4925G>A table: the reading beside step 46's count.
/// `kiv2_copies` is the count as its table prints it (`.` without one);
/// `why_not` is the reason the reads cannot be trusted (the aligner left
/// out the reads of MAPQ 0), which makes it not assessed whatever the counts.
pub fn null_allele(
    counts: Option<Counts>,
    kiv2_copies: &str,
    why_not: Option<&str>,
) -> Vec<(String, String)> {
    let m = &LPA_KIV2_4925;
    let a = assess(m, if why_not.is_some() { None } else { counts });
    let total = pycompat::num::parse_float(kiv2_copies)
        .ok()
        .map(|n| n + 2.0);
    let copies = copies_of(&a, total);
    let name = format!(
        "{} (a splice variant that lowers the Lp(a) made from its chromosome, Coassin 2017)",
        m.name
    );
    let st = a.status.as_str();
    let msg = if a.status == Status::NotAssessed {
        format!(
            "{name}: not assessed, {}",
            why_not.unwrap_or("the reads at its seven GRCh38 copies could not be read.")
        )
    } else {
        let mut s = format!(
            "{name}: {st}, {} of {} reads at its seven GRCh38 copies carry it",
            a.alt_reads, a.depth
        );
        if copies != "." {
            let t = total.expect("copies need the count");
            s += &format!(" (about {copies} of the {} copies counted)", fixed(t, 0));
        }
        if a.status == Status::TooLittleDepth {
            s += &format!(", under the {} reads a reading needs", m.min_depth);
        }
        s.push('.');
        match a.status {
            Status::Present => {
                s += " Carriers had a lower Lp(a) than non-carriers with the same number of repeats, most of all with few repeats. ";
                s += NULL_LIMIT;
            }
            Status::NotSeen => {
                s += " Not seen in the reads is not proof of absence: one carrying repeat among many gives few reads, and a short-read fraction never says which chromosome or repeat holds a variant."
            }
            _ => {
                s.push(' ');
                s += NULL_LIMIT;
            }
        }
        s
    };
    rows(&[
        ("marker", m.name.to_owned()),
        ("rsid", format!("{} (first KIV-2 copy)", m.dbsnp.0)),
        ("status", st.to_owned()),
        ("alt_reads", a.alt_reads),
        ("depth", a.depth),
        ("fraction", a.fraction),
        ("copies", copies),
        ("kiv2_copies", kiv2_copies.to_owned()),
        ("message", msg),
    ])
}

/// The reason a KIV-2 4925G>A reading cannot be trusted when reads of
/// MAPQ 0 make up less than `floor` of the KIV-2 depth: the reads that fit
/// several copies are the ones a carrying repeat gives.
pub fn mapq0_why_not(share: f64, floor: f64) -> Option<String> {
    (share < floor).then(|| {
        format!(
            "reads of mapping quality 0 make up {}% of the KIV-2 depth, so the aligner left out the reads that fit several copies, the ones a carrying repeat gives.",
            fixed(100.0 * share, 1)
        )
    })
}

/// SMN1's marker never proves 2+0, and its absence never rules it out.
pub const SMN1_NOT_PROOF: &str =
    "This marker does not prove 2+0, and its absence does not rule 2+0 out.";

/// Step 35's SMN1 marker table: the marker's reading beside the SMN1 psCN
/// of the SMN table (`smn`, read by [`crate::smn::read`]).
pub fn smn1_marker(counts: Option<Counts>, smn: &crate::smn::SmnRead) -> Vec<(String, String)> {
    let m = &SMN1_27134;
    let a = assess(m, counts);
    let copies = copies_of(&a, smn.agcn.map(|v| v as f64));
    let st = a.status.as_str();
    let head = if a.status == Status::NotAssessed {
        format!(
            "{} ({}): not assessed, the reads at the marker could not be read.",
            m.name, m.dbsnp.0
        )
    } else {
        let mut s = format!(
            "{} ({}): {st}, {} of {} reads carry G at the SMN1 and SMN2 copies of the site together",
            m.name, m.dbsnp.0, a.alt_reads, a.depth
        );
        if copies != "." {
            s += &format!(
                " (about {copies} of the {} SMN1 and SMN2 copies)",
                smn.agcn.expect("copies need the agCN")
            );
        }
        if a.status == Status::TooLittleDepth {
            s += &format!(", under the {} reads a reading needs", m.min_depth);
        }
        s + ". The reads cannot tell whether a G sits on SMN1 or SMN2."
    };
    let body = match smn.pscn {
        Some(2) => match a.status {
            Status::Present => " With two SMN1 copies, this marker goes with a higher chance that both sit on one chromosome (2+0, a silent carrier): in Luo 2014 it raised the residual carrier risk of a two-copy result (in African Americans from about 1 in 121 to 1 in 34); the figures depend on ancestry.".to_owned(),
            Status::NotSeen => " With two SMN1 copies and no marker, the chance of 2+0 is lower than for a two-copy result alone (in Luo 2014, in African Americans, about 1 in 396 against 1 in 121), not zero; the figures depend on ancestry.".to_owned(),
            _ => " With two SMN1 copies the marker would refine the chance of 2+0; here it gives no reading.".to_owned(),
        },
        Some(ps) => format!(
            " SMN1 has {ps} cop{} here; the marker refines the carrier risk of a two-copy result only.",
            if ps == 1 { "y" } else { "ies" }
        ),
        None => format!(
            " SMN1 copy number: {}; the marker refines the carrier risk of a two-copy result only.",
            smn.text
        ),
    };
    rows(&[
        ("marker", m.name.to_owned()),
        ("rsid", m.dbsnp.0.to_owned()),
        ("status", st.to_owned()),
        ("alt_reads", a.alt_reads),
        ("depth", a.depth),
        ("fraction", a.fraction),
        ("copies", copies),
        ("smn1_pscn", smn.text.clone()),
        ("agcn", crate::d(smn.agcn)),
        ("message", format!("{head}{body} {SMN1_NOT_PROOF}")),
    ])
}

fn rows(r: &[(&str, String)]) -> Vec<(String, String)> {
    r.iter()
        .map(|(k, v)| ((*k).to_owned(), v.clone()))
        .collect()
}

/// Problems of one marker (empty when it is consistent): the dbSNP record
/// is one of its copies with its ref and alt, on the plus strand and not
/// swapped; the context's middle base is the ref; no copy is listed twice.
pub fn marker_errors(m: &Marker) -> Vec<String> {
    let (rs, chrom, pos, r, a) = m.dbsnp;
    let comp = |b: &str| match b {
        "A" => "T",
        "C" => "G",
        "G" => "C",
        "T" => "A",
        _ => "",
    };
    let (mr, ma) = (
        (m.reference as char).to_string(),
        (m.alt as char).to_string(),
    );
    let mut errs = Vec::new();
    if chrom != m.chrom || !m.copies.contains(&pos) {
        errs.push(format!(
            "{}: dbSNP places {rs} at {chrom}:{pos}, not at one of the marker's copies",
            m.key
        ));
    }
    if (mr.as_str(), ma.as_str()) == (a, r) {
        errs.push(format!(
            "{}: ref and alt are swapped (marker {mr}>{ma}, dbSNP {r}>{a})",
            m.key
        ));
    } else if (mr.as_str(), ma.as_str()) == (comp(r), comp(a)) {
        errs.push(format!(
            "{}: the alleles are the minus-strand bases (marker {mr}>{ma}, dbSNP {r}>{a}): write the plus-strand bases",
            m.key
        ));
    } else if (mr.as_str(), ma.as_str()) != (r, a) {
        errs.push(format!("{}: marker {mr}>{ma}, dbSNP {r}>{a}", m.key));
    }
    let ctx = m.context.as_bytes();
    if ctx.len() % 2 != 1 || ctx[ctx.len() / 2] != m.reference {
        errs.push(format!(
            "{}: the middle base of the context is not the ref {mr}",
            m.key
        ));
    }
    let mut seen = std::collections::HashSet::new();
    if !m.copies.iter().all(|c| seen.insert(c)) {
        errs.push(format!("{}: a copy is listed twice", m.key));
    }
    errs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn markers_are_consistent() {
        for m in &MARKERS {
            assert!(marker_errors(m).is_empty(), "{:?}", marker_errors(m));
        }
        // A planted swap, a minus-strand pair and a copy listed twice are caught.
        let mut bad = SMN1_27134;
        bad.reference = b'G';
        bad.alt = b'T';
        assert!(marker_errors(&bad)[0].contains("swapped"));
        bad.reference = b'A';
        bad.alt = b'C';
        assert!(
            marker_errors(&bad)
                .iter()
                .any(|e| e.contains("minus-strand"))
        );
        bad = LPA_KIV2_4925;
        bad.copies = &[160_640_667, 160_640_667];
        assert!(
            marker_errors(&bad)
                .iter()
                .any(|e| e.contains("listed twice"))
        );
    }

    #[test]
    fn statuses() {
        let m = &LPA_KIV2_4925;
        assert_eq!(assess(m, Some((95, 5, 0))).status, Status::Present);
        assert_eq!(assess(m, Some((98, 2, 0))).status, Status::NotSeen);
        assert_eq!(assess(m, Some((996, 4, 0))).status, Status::Uncertain);
        assert_eq!(assess(m, Some((50, 5, 0))).status, Status::TooLittleDepth);
        assert_eq!(assess(m, None).status, Status::NotAssessed);
        assert_eq!(assess(m, Some((95, 5, 0))).fraction, "0.0500");
    }

    #[test]
    fn null_allele_copies_and_why_not() {
        let t = null_allele(Some((190, 10, 0)), "20.0", None);
        assert_eq!(t[6], ("copies".into(), "1.1".into()));
        assert!(
            t[8].1.contains("(about 1.1 of the 22 copies counted)"),
            "{}",
            t[8].1
        );
        let why = mapq0_why_not(0.1996, 0.2).unwrap();
        assert!(
            why.starts_with("reads of mapping quality 0 make up 20.0%"),
            "{why}"
        );
        let t = null_allele(Some((190, 10, 0)), "20.0", Some(&why));
        assert_eq!(t[2].1, "not assessed");
        assert!(mapq0_why_not(0.2, 0.2).is_none());
    }
}
