//! ClinVar's conflicting records that some labs call pathogenic: setup's
//! subset of them and step 06's table of the ones a sample carries.
//!
//! A record ClinVar labels `Conflicting_classifications_of_pathogenicity`
//! (`Conflicting_interpretations_of_pathogenicity` before 2024) is not in
//! setup's pathogenic subset, whose rule is `CLNSIG` holding `Pathogenic` or
//! `Likely_pathogenic`. When at least one lab submitted it as pathogenic or
//! likely pathogenic, a clinician would still want to see it. `CLNSIGCONF`
//! counts the submissions of each class, as in
//! `Pathogenic(2)|Uncertain_significance(1)`.
//!
//! - [`keep`] (`og clinvar conflicting-subset`): the records of the subset,
//!   `clinvar_conflicting_plp_chr.vcf.gz`. A record of the pathogenic subset
//!   is never kept, so the two subsets share no record.
//! - [`table`] (`og clinvar conflicting-hits`): one row per carried record,
//!   with the submissions, the stars and the condition, or a note that says
//!   the subset was not checked.

use std::io::{self, BufRead, Write};

use crate::text::{info_get, is_digits, nth, split};
use crate::{CallFloor, ReviewStatus, zygosity};

/// The columns of step 06's `<sample>_clinvar_conflicting.tsv`.
pub const COLUMNS: &str = "chrom\tpos\tref\talt\tgenotype\tgene\tclinvar_id\tstars\treview_status\tsubmissions\tplp_submissions\tcondition\tcondition_ids\tdp\tgq\tlow_confidence";

/// True when `CLNSIG` is ClinVar's conflicting label (one of the `|` or
/// `/` parts, under either name).
pub fn conflicting(clnsig: &[u8]) -> bool {
    clnsig.split(|&c| c == b'|' || c == b'/').any(|t| {
        t == b"Conflicting_classifications_of_pathogenicity"
            || t == b"Conflicting_interpretations_of_pathogenicity"
    })
}

/// True when setup's pathogenic subset keeps the record: its rule is
/// `CLNSIG~"Pathogenic" || CLNSIG~"Likely_pathogenic"` (a case-sensitive
/// substring), so `pathogenicity` of the conflicting label is no match.
pub fn in_plp_subset(clnsig: &[u8]) -> bool {
    memchr::memmem::find(clnsig, b"Pathogenic").is_some()
        || memchr::memmem::find(clnsig, b"Likely_pathogenic").is_some()
}

/// True for a pathogenic or likely pathogenic class of `CLNSIGCONF`:
/// `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`, each
/// also with `,_low_penetrance`, as setup's pathogenic subset counts them.
pub fn plp_class(class: &[u8]) -> bool {
    class.split(|&c| c == b'/').any(|p| {
        let p = p.strip_suffix(b",_low_penetrance").unwrap_or(p);
        p == b"Pathogenic" || p == b"Likely_pathogenic"
    })
}

/// The submissions `CLNSIGCONF` counts as pathogenic or likely pathogenic.
/// A class without a count, `Pathogenic` alone, is one submission.
pub fn plp_submissions(clnsigconf: &[u8]) -> u64 {
    let _ = clnsigconf;
    0
}

/// True when the subset keeps the record whose INFO column this is: a
/// conflicting record with at least one pathogenic or likely pathogenic
/// submission, that the pathogenic subset does not hold.
pub fn keep(info: &[u8]) -> bool {
    let _ = info;
    false
}

/// `og clinvar conflicting-subset`: the header of a ClinVar VCF as it is,
/// then the records [`keep`] keeps.
pub fn subset(mut input: impl BufRead, out: &mut impl Write) -> io::Result<()> {
    let mut line = Vec::new();
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Ok(());
        }
        let body = line.strip_suffix(b"\n").unwrap_or(&line);
        let body = body.strip_suffix(b"\r").unwrap_or(body);
        if body.is_empty() {
            continue;
        }
        if body.starts_with(b"#") || keep(nth(&split(body, b'\t'), 8)) {
            out.write_all(body)?;
            out.write_all(b"\n")?;
        }
    }
}

/// `_` as spaces, or `.` for an empty value.
fn words(v: &[u8]) -> Vec<u8> {
    if v.is_empty() || v == b"." {
        return b".".to_vec();
    }
    v.iter()
        .map(|&c| if c == b'_' { b' ' } else { c })
        .collect()
}

/// ClinVar's placeholders for a submission that names no condition.
fn placeholder(name: &[u8]) -> bool {
    let n = name.to_ascii_lowercase();
    matches!(n.as_slice(), b"" | b"not provided" | b"not specified")
}

/// The condition columns of a record: the names `CLNDN` gives (`_` read as
/// a space), joined by ` / `, ClinVar's placeholders left out; and the
/// database ids `CLNDISDB` gives for those conditions and for a condition
/// with ids and no name, comma-joined without repeats. `.` for none.
/// The run summary words a P/LP hit's condition the same way.
pub fn condition(clndn: &[u8], clndisdb: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let list = |v: &[u8]| -> Vec<Vec<u8>> {
        if v.is_empty() || v == b"." {
            Vec::new()
        } else {
            split(v, b'|').into_iter().map(<[u8]>::to_vec).collect()
        }
    };
    let (names, dbs) = (list(clndn), list(clndisdb));
    let mut named: Vec<Vec<u8>> = Vec::new();
    let mut ids: Vec<Vec<u8>> = Vec::new();
    for i in 0..names.len().max(dbs.len()) {
        let name: Vec<u8> = names
            .get(i)
            .map(|n| {
                n.iter()
                    .map(|&c| if c == b'_' { b' ' } else { c })
                    .collect::<Vec<u8>>()
                    .trim_ascii()
                    .to_vec()
            })
            .filter(|n| n != b".")
            .unwrap_or_default();
        if !name.is_empty() && placeholder(&name) {
            continue;
        }
        if !name.is_empty() {
            named.push(name);
        }
        for id in dbs.get(i).map(|d| split(d, b',')).unwrap_or_default() {
            let id = id.trim_ascii();
            if !id.is_empty() && id != b"." && !ids.iter().any(|x| x == id) {
                ids.push(id.to_vec());
            }
        }
    }
    let join = |v: Vec<Vec<u8>>, sep: &[u8]| {
        if v.is_empty() {
            b".".to_vec()
        } else {
            v.join(sep)
        }
    };
    (join(named, b" / "), join(ids, b","))
}

/// The row of one carried record of the annotated subset, under `floor`.
pub fn row(line: &[u8], floor: &CallFloor) -> Vec<u8> {
    let f = split(line, b'\t');
    let info = nth(&f, 8);
    let get = |k: &[u8]| info_get(info, k);
    let genes: Vec<&[u8]> = split(get(b"GENEINFO"), b'|')
        .into_iter()
        .filter(|g| !g.is_empty())
        .map(|g| nth(&split(g, b':'), 1))
        .collect();
    let gene = genes.join(&b',');
    let gene = if gene.is_empty() { b".".to_vec() } else { gene };
    let id = nth(&f, 3);
    let rev = get(b"CLNREVSTAT");
    let conf = get(b"CLNSIGCONF");
    let submissions = if conf.is_empty() {
        b".".to_vec()
    } else {
        split(conf, b'|').join(&b", "[..])
    };
    let (names, ids) = condition(get(b"CLNDN"), get(b"CLNDISDB"));
    let sample = split(nth(&f, 10), b':');
    let gt = nth(&sample, 1);
    let (mut dp, mut gq): (&[u8], &[u8]) = (b".", b".");
    for (i, k) in split(nth(&f, 9), b':').into_iter().enumerate() {
        let v = nth(&sample, i + 1);
        if k == b"DP" && is_digits(v) {
            dp = v;
        } else if k == b"GQ" && is_digits(v) {
            gq = v;
        }
    }
    let low = floor.low((gq != b".").then_some(gq), (dp != b".").then_some(dp));
    [
        nth(&f, 1),
        nth(&f, 2),
        nth(&f, 4),
        nth(&f, 5),
        &zygosity(gt),
        &gene,
        if id.is_empty() { &b"."[..] } else { id },
        ReviewStatus::parse(rev).stars().to_string().as_bytes(),
        &words(rev),
        &submissions,
        plp_submissions(conf).to_string().as_bytes(),
        &names,
        &ids,
        dp,
        gq,
        if low { b"yes" } else { b"no" },
    ]
    .join(&b'\t')
}

/// `og clinvar conflicting-hits`: the header, then one row per record of
/// the annotated subset hits. `skip` (why the subset was not checked) gives
/// a first line `# not checked: <why>`, the header and no row, so a missing
/// subset never reads as zero hits.
pub fn table(
    input: Option<impl BufRead>,
    skip: Option<&str>,
    floor: &CallFloor,
    out: &mut impl Write,
) -> io::Result<()> {
    if let Some(why) = skip {
        writeln!(out, "# not checked: {why}")?;
    }
    writeln!(out, "{COLUMNS}")?;
    let Some(mut input) = input else {
        return Ok(());
    };
    let mut line = Vec::new();
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Ok(());
        }
        let body = line.strip_suffix(b"\n").unwrap_or(&line);
        if !body.is_empty() && !body.starts_with(b"#") {
            out.write_all(&row(body, floor))?;
            out.write_all(b"\n")?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "##fileformat=VCFv4.1\n##fileDate=20261101\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn rec(pos: u32, info: &str) -> String {
        format!("chr22\t{pos}\t{pos}\tC\tT\t.\t.\t{info}\n")
    }

    fn run_subset(vcf: &str) -> String {
        let mut out = Vec::new();
        subset(vcf.as_bytes(), &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// The plan's fixture: a conflicting record two labs call pathogenic
    /// and one uncertain is kept; one with only benign and uncertain
    /// submissions is dropped; a plain pathogenic record stays out (the
    /// pathogenic subset holds it).
    #[test]
    fn the_subset_keeps_conflicting_records_with_a_plp_submission() {
        let kept = rec(
            100,
            "CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Pathogenic(2)|Uncertain_significance(1);CLNREVSTAT=criteria_provided,_conflicting_classifications",
        );
        let benign = rec(
            200,
            "CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Benign(1)|Uncertain_significance(3)",
        );
        let plain = rec(
            300,
            "CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter",
        );
        let got = run_subset(&format!("{HEAD}{kept}{benign}{plain}"));
        assert_eq!(got, format!("{HEAD}{kept}"));
    }

    #[test]
    fn the_older_label_and_every_plp_class_count() {
        for (conf, n) in [
            ("Likely_pathogenic(1)|Benign(4)", 1),
            ("Pathogenic/Likely_pathogenic(3)|Uncertain_significance(1)", 3),
            ("Pathogenic,_low_penetrance(1)|Likely_benign(1)", 1),
            ("Pathogenic(2)|Likely_pathogenic(1)|Benign(1)", 3),
            ("Likely_benign(2)|Uncertain_significance(1)", 0),
            ("Pathogenic|Benign(1)", 1),
            ("", 0),
        ] {
            assert_eq!(plp_submissions(conf.as_bytes()), n, "{conf}");
        }
        let old = "CLNSIG=Conflicting_interpretations_of_pathogenicity;CLNSIGCONF=Likely_pathogenic(1)|Benign(1)";
        assert!(keep(old.as_bytes()));
        // A conflicting record with no CLNSIGCONF cannot show a P/LP lab.
        assert!(!keep(
            b"CLNSIG=Conflicting_classifications_of_pathogenicity"
        ));
        // A CLNSIG the pathogenic subset also matches is left to it.
        assert!(!keep(
            b"CLNSIG=Conflicting_classifications_of_pathogenicity|Pathogenic;CLNSIGCONF=Pathogenic(1)"
        ));
        // Uncertain alone is not conflicting.
        assert!(!keep(
            b"CLNSIG=Uncertain_significance;CLNSIGCONF=Pathogenic(1)"
        ));
    }

    #[test]
    fn a_carried_record_in_sixteen_columns() {
        let line = "chr22\t100\t9900101\tC\tT\t.\tPASS\tGENEINFO=NF2:4771;CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Pathogenic(2)|Uncertain_significance(1);CLNREVSTAT=criteria_provided,_conflicting_classifications;CLNDN=Neurofibromatosis,_type_2|not_provided;CLNDISDB=MONDO:MONDO:0007039,OMIM:101000|MedGen:CN517202\tGT:GQ:DP\t0/1:20:30";
        let got = String::from_utf8(row(line.as_bytes(), &CallFloor::default())).unwrap();
        assert_eq!(
            got,
            "chr22\t100\tC\tT\thet\tNF2\t9900101\t1\tcriteria provided, conflicting classifications\tPathogenic(2), Uncertain_significance(1)\t2\tNeurofibromatosis, type 2\tMONDO:MONDO:0007039,OMIM:101000\t30\t20\tyes"
        );
        assert_eq!(got.split('\t').count(), COLUMNS.split('\t').count());
    }

    #[test]
    fn a_missing_subset_is_not_zero_rows() {
        let mut out = Vec::new();
        table(None::<&[u8]>, Some("run setup.sh --refresh clinvar"), &CallFloor::default(), &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            format!("# not checked: run setup.sh --refresh clinvar\n{COLUMNS}\n")
        );
        let mut out = Vec::new();
        table(Some(HEAD.as_bytes()), None, &CallFloor::default(), &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), format!("{COLUMNS}\n"));
    }

    #[test]
    fn conditions_drop_placeholders_and_repeats() {
        assert_eq!(
            condition(b"not_provided|Long_QT_syndrome_1", b"MedGen:CN517202|OMIM:192500,OMIM:192500"),
            (b"Long QT syndrome 1".to_vec(), b"OMIM:192500".to_vec())
        );
        assert_eq!(condition(b"", b""), (b".".to_vec(), b".".to_vec()));
        assert_eq!(
            condition(b"", b"OMIM:1"),
            (b".".to_vec(), b"OMIM:1".to_vec())
        );
    }
}
