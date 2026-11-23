//! Structural variants: Manta's and Delly's record counts, CNVpytor's
//! calls with the significance and quality filters, step 22's consensus
//! count, and step 05's class 4 and 5 calls with their support from the
//! consensus.

use super::roh::masked_overlap;
use super::util::{Obj, at, base, cells, count, int, lines, number, s, strs, table, vcf_lines};
use super::{Ctx, Found, Out};
use og_steps::sv::postfilter::{ConsensusSupport, Membership};
use pycompat::csv::DictRow;
use pycompat::json::Value;
use std::collections::HashMap;

/// A record count and the records that PASS.
fn count_vcf(path: &str) -> Result<(usize, usize), String> {
    let (mut total, mut pass) = (0, 0);
    vcf_lines(path, |l| {
        total += 1;
        if l.split(|&c| c == b'\t').nth(6) == Some(b"PASS") {
            pass += 1;
        }
    })?;
    Ok((total, pass))
}

fn totals(p: Found) -> Out {
    let (t, ps) = count_vcf(&p.path)?;
    Ok(Some((
        p,
        Obj::new()
            .set("total", count(t))
            .set("pass", count(ps))
            .val(),
    )))
}

/// Step 04 (Manta).
pub fn manta(c: &Ctx<'_>) -> Out {
    match c.find(&[
        "manta/results/variants/diploidSV.vcf.gz".to_owned(),
        "manta/*diploidSV.vcf.gz".to_owned(),
        "manta/**/diploidSV.vcf.gz".to_owned(),
    ]) {
        Some(p) => totals(p),
        None => Ok(None),
    }
}

/// Step 19 (Delly).
pub fn delly(c: &Ctx<'_>) -> Out {
    match c.find1(format!("delly/{}_sv.vcf.gz", c.s())) {
        Some(p) => totals(p),
        None => Ok(None),
    }
}

/// Step 22's consensus count.
pub fn sv_consensus(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("sv_merged/{}_sv_consensus.vcf.gz", c.s())) else {
        return Ok(None);
    };
    let (t, _) = count_vcf(&p.path)?;
    Ok(Some((p, Obj::new().set("consensus", count(t)).val())))
}

/// A CNVpytor call is kept in the filtered count when fewer than half of
/// its reads have mapping quality 0 and fewer than half of its bases are N,
/// and it touches no masked band.
const CNV_MAX_Q0: f64 = 0.5;
const CNV_MAX_PN: f64 = 0.5;

/// `(\S+):(\d+)-(\d+)` over a whole cell.
fn span(x: &str) -> Option<(&str, i128, i128)> {
    let (c, rest) = x.rsplit_once(':')?;
    let (a, b) = rest.split_once('-')?;
    let digits = |v: &str| !v.is_empty() && v.bytes().all(|d| d.is_ascii_digit());
    if c.is_empty() || c.chars().any(char::is_whitespace) || !digits(a) || !digits(b) {
        return None;
    }
    Some((c, a.parse().ok()?, b.parse().ok()?))
}

/// Step 18 (CNVpytor): the calls, the significant ones (e-value under
/// 0.01), and those kept after the quality and band filters.
pub fn cnvpytor(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("cnvpytor/{}_cnvs.txt", c.s())) else {
        return Ok(None);
    };
    let (mut total, mut dels, mut dups, mut sig, mut kept, mut kept_y) = (0, 0, 0, 0, 0, 0);
    let (mut q0s, mut gaps, mut bands) = (0, 0, 0);
    for line in lines(&p.path)? {
        let c = pycompat::text::split(&line, None, None);
        if c.is_empty() {
            continue;
        }
        total += 1;
        dels += usize::from(c[0] == "deletion");
        dups += usize::from(c[0] == "duplication");
        let m = if c.len() > 4 { span(c[1]) } else { None };
        match c.get(4).map(|v| pycompat::num::parse_float(v)) {
            Some(Ok(e)) if e.partial_cmp(&0.01).is_none_or(|o| o.is_lt()) => {}
            _ => continue,
        }
        sig += 1;
        let (q0, pn) = if c.len() > 9 {
            (number(c[8]), number(c[9]))
        } else {
            (None, None)
        };
        if q0.is_some_and(|q| q >= CNV_MAX_Q0) {
            q0s += 1;
        } else if pn.is_some_and(|q| q >= CNV_MAX_PN) {
            gaps += 1;
        } else if let Some((ch, a, b)) = m
            && masked_overlap(ch, a, b) > 0
        {
            bands += 1;
        } else {
            kept += 1;
            kept_y += usize::from(m.is_some_and(|(ch, ..)| ch == "chrY" || ch == "Y"));
        }
    }
    let data = Obj::new()
        .set("total", count(total))
        .set("deletions", count(dels))
        .set("duplications", count(dups))
        .set("significant", count(sig))
        .set("significant_filtered", count(kept))
        .set("filtered_chry", count(kept_y))
        .set(
            "dropped",
            Obj::new()
                .set("q0", count(q0s))
                .set("gap", count(gaps))
                .set("band", count(bands))
                .val(),
        );
    Ok(Some((p, data.val())))
}

/// Takes the chrY calls out of the filtered count of a declared female;
/// with no declared sex they stay, and the data says so.
pub fn cnvpytor_sex_filter(v: &mut Value, declared_sex: Option<&str>) {
    let Value::Object(m) = v else { return };
    let (Some(Value::Int(kept)), Some(Value::Int(y))) = (
        m.get("significant_filtered").cloned(),
        m.get("filtered_chry").cloned(),
    ) else {
        return;
    };
    let declared = declared_sex.filter(|x| !x.is_empty());
    m.insert("declared_sex".into(), declared.map_or(Value::Null, s));
    if declared == Some("female") {
        if let Some(Value::Object(d)) = m.get_mut("dropped") {
            d.insert("chry_female".into(), Value::Int(y));
        }
        m.insert("significant_filtered".into(), Value::Int(kept - y));
    }
}

// --- step 05, AnnotSV class 4 and 5 ---------------------------------------------

/// Support and callers from the SURVIVOR event containing the source call.
/// A row in no event has its own caller alone. Both values are null when
/// the membership is not known ([`Membership::Unknown`]) and without a
/// consensus file.
fn consensus_support(rows: &[DictRow], path: Option<&str>) -> Result<Vec<(Value, Value)>, String> {
    let Some(path) = path.filter(|p| crate::is_file(p)) else {
        return Ok(rows.iter().map(|_| (Value::Null, Value::Null)).collect());
    };
    let mut consensus = ConsensusSupport::default();
    for line in lines(path)? {
        consensus.add_line(&line);
    }
    Ok(rows
        .iter()
        .map(|row| {
            let cell = |key| row.get(key).flatten().unwrap_or_default();
            match consensus.membership(
                cell("SV_chrom"),
                cell("SV_type"),
                cell("Caller"),
                cell("ID"),
            ) {
                Membership::Member(callers) => (count(callers.len()), s(&callers.join(","))),
                Membership::NotMember => (int(1), s(&super::util::or(row, "Caller", "."))),
                Membership::Unknown => (Value::Null, Value::Null),
            }
        })
        .collect())
}

/// `{FORMAT key: value}` of the sample column (the column after FORMAT,
/// unless it is `Annotation_mode`).
fn sample_format(r: &DictRow, head: &[String]) -> HashMap<String, String> {
    let keys: Vec<&str> = r.iter().map(|(k, _)| k).collect();
    let _ = head;
    let Some(i) = keys.iter().position(|k| *k == "FORMAT").map(|i| i + 1) else {
        return HashMap::new();
    };
    let Some(col) = keys.get(i).filter(|k| **k != "Annotation_mode") else {
        return HashMap::new();
    };
    let fmt = r.get("FORMAT").flatten().unwrap_or_default();
    let val = r.get(col).flatten().unwrap_or_default();
    fmt.split(':')
        .zip(val.split(':'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect()
}

/// Step 05's class 4 and 5 SVs, each with the callers that support it from
/// step 22's consensus; a call step 05's rules moved to class 3 listed
/// apart.
pub fn annotsv(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("annotsv/{}_sv_class4_5_checked.tsv", c.s())) else {
        return Ok(None);
    };
    let size = std::fs::metadata(&p.path).map_or(0, |m| m.len());
    let t = if size > 0 {
        table(&p.path, false)?
    } else {
        super::util::Table {
            head: Vec::new(),
            rows: Vec::new(),
        }
    };
    let need = [
        "SV_chrom",
        "SV_start",
        "SV_end",
        "SV_type",
        "ACMG_class",
        "Caller",
        "Exon_edge",
        "Population_AF",
        "Common",
    ];
    if !t.rows.is_empty() && !need.iter().all(|k| t.has(k)) {
        let mut sorted = need.to_vec();
        sorted.sort_unstable();
        return Err(format!(
            "{} lacks a column of [{}]",
            base(&p.path),
            sorted
                .iter()
                .map(|k| crate::repr(k))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    let mut cons = c.find1(format!("sv_merged/{}_sv_consensus.vcf.gz", c.s()));
    let mut cons_note = Value::Null;
    let support = match consensus_support(&t.rows, cons.as_ref().map(|f| f.path.as_str())) {
        Ok(x) => x,
        Err(e) => {
            let f = cons.take().map(|f| f.path).unwrap_or_default();
            let why = e.split_once(": ").map_or(e.as_str(), |x| x.1);
            cons_note = s(&format!("{} could not be read: {why}", base(&f)));
            consensus_support(&t.rows, None)?
        }
    };
    let mut keep = Vec::new();
    let mut moved = Vec::new();
    for (r, (sp, callers)) in t.rows.iter().zip(support) {
        // A row cut short reads None there, as Python's DictReader gives it.
        let cell =
            |k: &str| -> Result<String, String> { Ok(at(r, k)?.unwrap_or("None").to_owned()) };
        let cellv = |k: &str| -> Result<Value, String> { Ok(super::util::opt_s(at(r, k)?)) };
        let af = match at(r, "Population_AF")? {
            None | Some("") => Value::Null,
            Some(v) => pycompat::num::parse_float(v).map_or(Value::Null, Value::Float),
        };
        let genes: Vec<String> = super::util::or(r, "Gene_name", "")
            .split(';')
            .filter(|g| !g.is_empty())
            .map(str::to_owned)
            .collect();
        let fmt = sample_format(r, &t.head);
        let fmt_get = |k: &str| {
            fmt.get(k)
                .filter(|v| !v.is_empty())
                .map_or(Value::Null, |v| s(v))
        };
        let acmg = cell("ACMG_class")?;
        let support_callers = match &callers {
            Value::Str(x) if !x.is_empty() => strs(x.split(',')),
            _ => Value::Array(Vec::new()),
        };
        let row = Obj::new()
            .set(
                "sv",
                s(&format!(
                    "{}:{}-{}",
                    cell("SV_chrom")?,
                    cell("SV_start")?,
                    cell("SV_end")?
                )),
            )
            .set("type", cellv("SV_type")?)
            .set("gt", fmt_get("GT"))
            .set("ps", fmt_get("PS"))
            .set("length", s(&super::util::or(r, "SV_length", "")))
            .set("acmg_class", cellv("ACMG_class")?)
            .set("genes", strs(&genes))
            .set("caller", cellv("Caller")?)
            .set("support", sp)
            .set("support_callers", support_callers)
            .set("exon_edge", cellv("Exon_edge")?)
            .set(
                "exon_edge_detail",
                s(&super::util::or(r, "Exon_edge_detail", "")),
            )
            .set("population_af", af)
            .set(
                "population_source",
                s(&super::util::or(r, "Population_source", "")),
            )
            .set("common", Value::Bool(at(r, "Common")? == Some("yes")));
        // Step 05 writes Benign_match when AnnotSV gives the benign records'
        // coordinates; a table without it gives no field.
        let row = if t.has("Benign_match") {
            row.set("benign_match", cellv("Benign_match")?)
        } else {
            row
        };
        let row = row
            .set(
                "annotsv_class",
                match r.get("AnnotSV_class").flatten().filter(|v| !v.is_empty()) {
                    Some(v) => s(v),
                    None => cellv("ACMG_class")?,
                },
            )
            .set("rule", s(&super::util::or(r, "Rule", "")))
            .set("rule_detail", s(&super::util::or(r, "Rule_detail", "")));
        if matches!(acmg.as_str(), "4" | "5") {
            keep.push(row.val());
        } else {
            moved.push(row.val());
        }
    }
    let field = |v: &Value, k: &str| -> Value {
        match v {
            Value::Object(m) => m.get(k).cloned().unwrap_or(Value::Null),
            _ => Value::Null,
        }
    };
    let single = keep
        .iter()
        .filter(|x| field(x, "support") == Value::Int(1))
        .count();
    let edge = keep
        .iter()
        .filter(|x| field(x, "exon_edge") == s("yes"))
        .count();
    let edge_checked = !keep
        .iter()
        .any(|x| field(x, "exon_edge") == s("not checked"));
    let common = keep
        .iter()
        .filter(|x| field(x, "common") == Value::Bool(true))
        .count();
    let data = Obj::new()
        .set("count", count(keep.len()))
        .set("rows", Value::Array(keep))
        .set("downgraded", Value::Array(moved))
        .set("consensus", cons.map_or(Value::Null, |f| s(&f.rel)))
        .set("consensus_note", cons_note)
        .set("single_caller", count(single))
        .set("exon_edge", count(edge))
        .set("exon_edge_checked", Value::Bool(edge_checked))
        .set("common", count(common));
    let _ = cells;
    Ok(Some((p, data.val())))
}

#[cfg(test)]
mod tests {
    use super::*;

    const REPORT_HEADER: &str = "ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tACMG_class\tCaller\tGene_name\tExon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\tAnnotSV_class\tRule\tRule_detail\n";
    const REPORT_ROWS: &str = concat!(
        "MantaDEL:1:0\t20\t1000\t5000\tDEL\t5\tmanta\tG1;G2;G3;G4;G5\tno\t\t\t\tno\t5\t\t\n",
        "NEARBY\t20\t3000\t7000\tDEL\t4\tmanta\tFIXG1\tyes\texon 1 and exon 3 edges\t0.05\tfixture\tyes\t4\t\t\n",
        "OTHER\t20\t3100\t7100\tDEL\t4\tmanta\tG9\tno\t\t\t\tno\t4\t\t\n",
    );
    const REPORT_LAST: &str = "RULE\t20\t50000\t55000\tDEL\t4\tmanta\tRULEG\tno\t\t\t\tno\t4\t\t\n";
    const REPORT_MOVED: &str = "RULE\t20\t50000\t55000\tDEL\t3\tmanta\tRULEG\tno\t\t\t\tno\t4\tig_tr_locus\tinside the IG locus\n";
    const REPORT_CONSENSUS: &str = concat!(
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\n",
        "chr20\t3000\tCLUSTER\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=7000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:MantaDEL_1_0\t0/1:DEL:D1\n",
    );

    fn report_field<'a>(value: &'a Value, key: &str) -> &'a Value {
        let Value::Object(fields) = value else {
            panic!("expected an object for {key}: {value:?}");
        };
        &fields[key]
    }

    fn report_rows<'a>(data: &'a Value, key: &str) -> &'a [Value] {
        let Value::Array(rows) = report_field(data, key) else {
            panic!("expected rows in {key}");
        };
        rows
    }

    fn report_counts(data: &Value) -> [&Value; 4] {
        ["count", "single_caller", "exon_edge", "common"].map(|key| report_field(data, key))
    }

    #[test]
    fn annotsv_report_payload_tracks_membership_and_rule_downgrades() {
        use crate::collect::{Group, Input, collect};
        use crate::testutil::Dir;

        let dir = Dir::new("annotsv-report-support");
        dir.put_gz("S/sv_merged/S_sv_consensus.vcf.gz", REPORT_CONSENSUS);
        let paired = |last: &str| {
            dir.put(
                "S/annotsv/S_sv_class4_5_checked.tsv",
                format!("{REPORT_HEADER}{REPORT_ROWS}{last}"),
            );
            let summary = collect(&Input {
                sample: "S".into(),
                sample_dir: dir.path("S"),
                declared_sex: None,
                assets: dir.path("assets"),
                group: Group::C,
                tables: dir.path("bin"),
            })
            .unwrap()
            .to_json();
            let section = report_field(report_field(&summary, "sections"), "annotsv");
            assert_eq!(report_field(section, "state"), &s("ok"));
            report_field(section, "data").clone()
        };
        let full = paired(REPORT_LAST);
        let moved = paired(REPORT_MOVED);
        let rows = report_rows(&full, "rows");
        let member = &rows[0];
        let nearby = &rows[1];
        let two_callers = report_field(member, "support") == &int(2)
            && report_field(member, "support_callers") == &strs(["manta", "delly"]);
        let one_caller = report_field(nearby, "support") == &int(1)
            && report_field(nearby, "support_callers") == &strs(["manta"]);
        let [total, single, edge, common] = report_counts(&full);
        let [after_total, after_single, after_edge, after_common] = report_counts(&moved);
        let downgraded = &report_rows(&moved, "downgraded")[0];

        // These checks read summary.json fields only; the text and HTML reports
        // render them, and og-render's own tests cover the labels. Keep every
        // migrated assertion visible in the control, even when several fail.
        let checks = [
            (
                "summary support and support_callers of a member",
                two_callers,
            ),
            (
                "summary single_caller, exon_edge and common counts",
                (single, edge, common) == (&int(3), &int(1), &int(1)),
            ),
            (
                "summary count with the single_caller, exon_edge and common counts",
                (total, single, edge, common) == (&int(4), &int(3), &int(1), &int(1)),
            ),
            (
                "summary sv, type, acmg_class and genes of the two-caller deletion",
                two_callers
                    && report_field(member, "sv") == &s("20:1000-5000")
                    && report_field(member, "type") == &s("DEL")
                    && report_field(member, "acmg_class") == &s("5")
                    && report_field(member, "genes") == &strs(["G1", "G2", "G3", "G4", "G5"]),
            ),
            (
                "summary support of a member and of a nearby non-member, with its caller",
                two_callers && one_caller && report_field(nearby, "caller") == &s("manta"),
            ),
            (
                "summary counts with exon_edge_checked",
                (total, single, edge, common) == (&int(4), &int(3), &int(1), &int(1))
                    && report_field(&full, "exon_edge_checked") == &Value::Bool(true),
            ),
            (
                "summary counts after excluding a downgraded call",
                (after_single, after_edge, after_common) == (&int(2), &int(1), &int(1)),
            ),
            (
                "summary counts and the downgraded row after a rule moved it",
                (after_total, after_single, after_edge, after_common)
                    == (&int(3), &int(2), &int(1), &int(1))
                    && report_field(downgraded, "sv") == &s("20:50000-55000")
                    && report_field(downgraded, "annotsv_class") == &s("4")
                    && report_field(downgraded, "acmg_class") == &s("3")
                    && report_field(downgraded, "rule") == &s("ig_tr_locus")
                    && report_rows(&moved, "rows")
                        .iter()
                        .all(|row| report_field(row, "sv") != &s("20:50000-55000")),
            ),
            (
                "summary exon_edge_detail, population and common fields of a row",
                report_field(nearby, "exon_edge_detail") == &s("exon 1 and exon 3 edges")
                    && report_field(nearby, "population_af") == &Value::Float(0.05)
                    && report_field(nearby, "population_source") == &s("fixture")
                    && report_field(nearby, "common") == &Value::Bool(true),
            ),
        ];
        let failed: Vec<&str> = checks
            .into_iter()
            .filter_map(|(name, passed)| (!passed).then_some(name))
            .collect();
        assert!(failed.is_empty(), "failed report checks: {failed:#?}");
    }

    #[test]
    fn the_benign_match_reaches_the_summary_only_when_step_05_wrote_it() {
        use crate::collect::{Group, Input, collect};
        use crate::testutil::Dir;

        let dir = Dir::new("annotsv-benign-match");
        let rows = |table: String| {
            dir.put("S/annotsv/S_sv_class4_5_checked.tsv", table);
            let summary = collect(&Input {
                sample: "S".into(),
                sample_dir: dir.path("S"),
                declared_sex: None,
                assets: dir.path("assets"),
                group: Group::C,
                tables: dir.path("bin"),
            })
            .unwrap()
            .to_json();
            let section = report_field(report_field(&summary, "sections"), "annotsv");
            report_rows(report_field(section, "data"), "rows").to_vec()
        };
        let head = REPORT_HEADER.replace("Common\t", "Common\tBenign_match\t");
        let with = rows(format!(
            "{head}{}{}",
            "NESTED\t20\t3000\t3500\tDEL\t4\tmanta\tG1\tno\t\t0.05\tgnomAD-SV_v4\tyes\tcontained\t4\t\t\n",
            "BLANK\t20\t9000\t9500\tDEL\t4\tmanta\tG2\tno\t\t\tno frequency data\tno\tnone\t4\t\t\n",
        ));
        let field = |r: &Value, k: &str| report_field(r, k).clone();
        assert_eq!(
            (field(&with[0], "benign_match"), field(&with[0], "common")),
            (s("contained"), Value::Bool(true)),
            "Common keeps its meaning beside the match"
        );
        assert_eq!(
            (
                field(&with[1], "population_af"),
                field(&with[1], "population_source"),
                field(&with[1], "benign_match")
            ),
            (Value::Null, s("no frequency data"), s("none"))
        );
        let without = rows(format!("{REPORT_HEADER}{REPORT_LAST}"));
        let Value::Object(m) = &without[0] else {
            panic!("a row is an object")
        };
        assert!(!m.contains_key("benign_match"), "{m:?}");
    }

    #[test]
    fn reports_and_two_hit_rules_use_the_original_member_id() {
        let d = std::env::temp_dir().join(format!("og-summary-sv-support-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let checked = d.join("checked.tsv").to_string_lossy().into_owned();
        let consensus = d.join("consensus.vcf").to_string_lossy().into_owned();
        std::fs::write(
            &checked,
            concat!(
                "ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tACMG_class\tCaller\tExon_edge\tPopulation_AF\tCommon\n",
                "MantaDEL:1:0\t20\t1000\t5000\tDEL\t4\tmanta\tno\t\tno\n",
                "UNLISTED\t20\t2000\t6000\tDEL\t4\tmanta\tno\t\tno\n",
                ".\t20\t2000\t6000\tDEL\t4\tmanta\tno\t\tno\n",
            ),
        )
        .unwrap();
        std::fs::write(
            &consensus,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\n",
                "chr20\t2000\tCLUSTER\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=6000;SUPP=2;SUPP_VEC=11\tGT:TY:ID\t0/1:DEL:MantaDEL_1_0\t0/1:DEL:D1\n",
            ),
        )
        .unwrap();
        let rows = table(&checked, false).unwrap().rows;
        let report = consensus_support(&rows, Some(&consensus)).unwrap();
        let two_hits = crate::collect::two_hits::read_annotsv(
            Some((checked, "checked.tsv".to_owned())),
            Some(consensus),
        )
        .unwrap()
        .unwrap();
        // A row with ID `.` cannot be matched to a member: its support is
        // not known (null), as without a consensus, never single caller.
        assert_eq!(
            (
                &report[..],
                two_hits[0].support,
                two_hits[1].support,
                two_hits[2].support
            ),
            (
                &[
                    (int(2), s("manta,delly")),
                    (int(1), s("manta")),
                    (Value::Null, Value::Null)
                ][..],
                Some(2),
                Some(1),
                None
            ),
        );
        std::fs::remove_dir_all(d).unwrap();
    }
}
