//! Step 42's traits section: the traits table, its genotype-rules lines,
//! and the APOL1 G2 deletion read from the sample's VCF beside the APOL1
//! rule, whose single-base sites read G1 only.

use super::util::{Obj, base, count, ragged, s, table};
use super::{Ctx, Out};
use pycompat::json::Value;

/// The columns of step 42's genotype-rules table.
pub const RULE_OUT_COLUMNS: [&str; 15] = [
    "rule",
    "gene",
    "trait",
    "rsids",
    "chrom",
    "pos",
    "genotypes",
    "gq",
    "phase",
    "result",
    "status",
    "reading",
    "ancestry",
    "note",
    "pmids",
];
/// The rules whose line the traits table holds already.
const TRAITS_TABLE_RULES: [&str; 1] = ["APOE"];
const TRAIT_COLUMNS: [&str; 9] = [
    "rsid", "gene", "trait", "genotype", "status", "reading", "ancestry", "note", "pmids",
];

pub fn traits(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("traits/{}_traits.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, true)?;
    if !t.first_has(&TRAIT_COLUMNS) {
        return Err(format!(
            "{} has no rows or lacks a column of {}",
            base(&p.path),
            TRAIT_COLUMNS.join(", ")
        ));
    }
    if let Some(i) = t.rows.iter().position(|r| ragged(r, &TRAIT_COLUMNS)) {
        return Err(format!(
            "{} line {} does not have one field per column",
            base(&p.path),
            i + 2
        ));
    }
    let cell = |r: &pycompat::csv::DictRow, k: &str| r.get(k).flatten().unwrap_or("").to_owned();
    let row = |r: &pycompat::csv::DictRow| {
        TRAIT_COLUMNS
            .iter()
            .fold(Obj::new(), |o, k| o.set(k, s(&cell(r, k))))
            .val()
    };
    let apoe = |r: &&pycompat::csv::DictRow| cell(r, "gene") == "APOE";
    let mut out = Obj::new()
        .set(
            "rows",
            Value::Array(t.rows.iter().filter(|r| !apoe(r)).map(row).collect()),
        )
        .set(
            "apoe",
            Value::Array(t.rows.iter().filter(apoe).map(row).collect()),
        );
    for st in ["called", "not_covered", "other_allele"] {
        out.put(
            st,
            count(
                t.rows
                    .iter()
                    .filter(|r| !apoe(r) && cell(r, "status") == st)
                    .count(),
            ),
        );
    }
    out.put("rules", genotype_rules(c)?);
    Ok(Some((p, out.val())))
}

/// Step 42's genotype-rules lines, `null` without the table; the APOL1 line
/// gets the G2 reading.
fn genotype_rules(c: &Ctx<'_>) -> Result<Value, String> {
    let Some(p) = c.find1(format!("traits/{}_genotype_rules.tsv", c.s())) else {
        return Ok(Value::Null);
    };
    let t = table(&p.path, true)?;
    if !t.rows.is_empty() && RULE_OUT_COLUMNS.iter().any(|k| !t.has(k)) {
        return Err(format!(
            "{} lacks a column of {}",
            base(&p.path),
            RULE_OUT_COLUMNS.join(", ")
        ));
    }
    if let Some(i) = t.rows.iter().position(|r| ragged(r, &RULE_OUT_COLUMNS)) {
        return Err(format!(
            "{} line {} does not have one field per column",
            base(&p.path),
            i + 2
        ));
    }
    let mut out = Vec::new();
    for r in &t.rows {
        let cell = |k: &str| r.get(k).flatten().unwrap_or("").to_owned();
        if TRAITS_TABLE_RULES.contains(&cell("rule").as_str()) {
            continue;
        }
        let mut o = RULE_OUT_COLUMNS
            .iter()
            .fold(Obj::new(), |o, k| o.set(k, s(&cell(k))));
        if cell("rule") == apol1::G2_RULE {
            let g2 = apol1::read(c);
            let both = apol1::combine(&cell("result"), &cell("status"), &g2);
            let mut v = g2.json();
            v.put("result", s(&both.0));
            v.put("reading", s(&both.1));
            o.put("g2", v.val());
        }
        out.push(o.val());
    }
    Ok(Value::Array(out))
}

/// APOL1 G2 (rs71785313), a 6-base deletion: read from the VCF records over
/// its span, and put together with the G1 result of step 42's APOL1 rule.
pub mod apol1 {
    use super::*;

    pub const G2_RULE: &str = "APOL1";
    const CHROM: &str = "chr22";
    const POS: i64 = 36_265_995;
    const REF: &str = "AATAATT";
    const ALT: &str = "A";
    /// GRCh38 chr22:36265975-36266025, the plus strand.
    const CONTEXT_START: i64 = 36_265_975;
    const CONTEXT: &str = "AGAAGCTAAACATTCTCAACAATAATTATAAGATTCTGCAGGCGGACCAAG";
    const CONTEXT_END: i64 = CONTEXT_START + CONTEXT.len() as i64 - 1;
    /// The anchor base of the left-aligned record to the end of the span a
    /// 6-base deletion there can be written over.
    pub const SPAN: (i64, i64) = (36_265_995, 36_266_005);
    const TWO_COPIES: &str = "the high-risk genotype of Genovese 2010, which went with focal segmental glomerulosclerosis and hypertension-attributed kidney failure in African Americans. A risk factor, not a diagnosis.";

    /// What `read_g2` gave, with the file it read.
    #[derive(Debug, Clone, PartialEq)]
    pub struct G2 {
        pub status: &'static str,
        pub copies: Option<i64>,
        pub gq: String,
        pub text: String,
        /// The VCF read, relative to the folder; `None` with no VCF. A
        /// read error writes it before the text.
        pub source: Option<String>,
        pub source_first: bool,
    }

    impl G2 {
        pub fn json(&self) -> Obj {
            let o = Obj::new()
                .set("status", s(self.status))
                .set(
                    "copies",
                    self.copies.map_or(Value::Null, |c| Value::Int(c.into())),
                )
                .set("gq", s(&self.gq));
            let src = super::super::util::opt_s(self.source.as_deref());
            if self.source_first {
                o.set("source", src).set("text", s(&self.text))
            } else {
                o.set("text", s(&self.text)).set("source", src)
            }
        }
    }

    /// The context with `r` at `pos` replaced by `a`, or `None` when `r` is
    /// not the GRCh38 bases there.
    pub fn apply(pos: i64, r: &str, a: &str) -> Option<String> {
        let i = pos - CONTEXT_START;
        if i < 0 || pos + r.len() as i64 - 1 > CONTEXT_END {
            return None;
        }
        let i = i as usize;
        if CONTEXT.get(i..i + r.len()) != Some(r) {
            return None;
        }
        Some(format!("{}{a}{}", &CONTEXT[..i], &CONTEXT[i + r.len()..]))
    }

    fn haplotype() -> Option<String> {
        apply(POS, REF, ALT)
    }

    fn digits(x: &str) -> bool {
        pycompat::text::isdigit(x)
    }

    fn record_end(r: &[String]) -> Option<i64> {
        if let Some(info) = r.get(7) {
            for (i, _) in info.match_indices("END=") {
                if i == 0 || info.as_bytes()[i - 1] == b';' {
                    let d: String = info[i + 4..]
                        .chars()
                        .take_while(char::is_ascii_digit)
                        .collect();
                    if !d.is_empty() {
                        return d.parse().ok();
                    }
                }
            }
        }
        Some(r[1].parse::<i64>().ok()? + r[3].chars().count() as i64 - 1)
    }

    fn overlaps(r: &[String]) -> bool {
        r.len() > 4
            && r[0] == CHROM
            && digits(&r[1])
            && r[1].parse::<i64>().is_ok_and(|p| p <= SPAN.1)
            && record_end(r).is_some_and(|e| e >= SPAN.0)
    }

    fn symbolic(a: &str) -> bool {
        a == "*" || a == "." || (a.len() >= 2 && a.starts_with('<') && a.ends_with('>'))
    }

    /// G2 from the records over [`SPAN`]; `gvcf` says a base with no
    /// record has no reads.
    pub fn read_g2(records: &[Vec<String>], gvcf: bool) -> G2 {
        let g2_hap = haplotype();
        let (mut copies, mut covered, mut gqs, mut problems) = (
            0i64,
            Vec::<(i64, i64)>::new(),
            Vec::<String>::new(),
            Vec::<String>::new(),
        );
        let (mut seen, mut g2_records) = (false, 0);
        'records: for r in records {
            if !overlaps(r) {
                continue;
            }
            if r.len() < 10 {
                problems.push("a record over G2 has no sample column".into());
                continue;
            }
            seen = true;
            let pos: i64 = r[1].parse().unwrap_or(0);
            let end = record_end(r).unwrap_or(pos);
            let rf = r[3].as_str();
            let alleles: Vec<&str> = std::iter::once(rf).chain(r[4].split(',')).collect();
            let mut fmt = indexmap::IndexMap::new();
            for (k, v) in r[8].split(':').zip(r[9].split(':')) {
                fmt.insert(k, v);
            }
            let field = |k: &str| fmt.get(k).copied().unwrap_or(".").to_owned();
            let gt = field("GT");
            let gt: Vec<&str> = gt.split(['/', '|']).collect();
            if gt.iter().any(|g| !digits(g)) {
                problems.push(format!("a no-call at {}:{pos}", r[0]));
                continue;
            }
            let idx: Vec<usize> = gt.iter().map(|g| g.parse().unwrap_or(usize::MAX)).collect();
            if idx.iter().any(|i| *i >= alleles.len()) {
                problems.push(format!(
                    "a genotype at {}:{pos} names an allele the record lacks",
                    r[0]
                ));
                continue;
            }
            if idx.iter().all(|i| *i == 0) {
                covered.push((pos, end));
                gqs.push(field("GQ"));
                continue;
            }
            if r[6] != "PASS" && r[6] != "." {
                problems.push(format!(
                    "a call at {}:{pos} that failed the caller's filter ({})",
                    r[0], r[6]
                ));
                continue;
            }
            let mut n = 0;
            let mut stopped = false;
            for &i in &idx {
                if i == 0 {
                    continue;
                }
                let a = alleles[i];
                let inside = pos >= CONTEXT_START && end <= CONTEXT_END;
                if inside && apply(pos, rf, rf).is_none() {
                    problems.push(format!(
                        "REF {rf} at {}:{pos} is not the GRCh38 base (is the VCF on GRCh38?)",
                        r[0]
                    ));
                    stopped = true;
                    break;
                }
                let hap = if symbolic(a) || !inside {
                    None
                } else {
                    apply(pos, rf, a)
                };
                if hap.is_some() && hap == g2_hap {
                    n += 1;
                } else {
                    problems.push(format!(
                        "an allele other than G2 inside its span ({}:{pos} {rf}>{a})",
                        r[0]
                    ));
                    stopped = true;
                    break;
                }
            }
            if stopped {
                continue;
            }
            if n > 0 {
                g2_records += 1;
                if g2_records > 1 {
                    problems.push("more than one record calls the G2 deletion".into());
                    break 'records;
                }
            }
            copies += n;
            covered.push((pos, end));
            gqs.push(field("GQ"));
        }
        let gq = if gqs.is_empty() {
            ".".to_owned()
        } else {
            gqs.join(",")
        };
        let g2 = |status, copies, gq: String, text: String| G2 {
            status,
            copies,
            gq,
            text,
            source: None,
            source_first: false,
        };
        if let Some(p) = problems.first() {
            return g2(
                "not_read",
                None,
                gq,
                format!("G2 (rs71785313) not read: {p}."),
            );
        }
        if copies > 0 {
            return g2(
                "called",
                Some(copies),
                gq,
                format!(
                    "G2 (rs71785313): {copies} cop{} of the deletion called.",
                    if copies > 1 { "ies" } else { "y" }
                ),
            );
        }
        if (SPAN.0..=SPAN.1).all(|p| covered.iter().any(|(a, b)| *a <= p && p <= *b)) {
            return g2(
                "called",
                Some(0),
                gq,
                "G2 (rs71785313): the reference over the whole deletion span, no G2 copy.".into(),
            );
        }
        if !gvcf && !seen {
            return g2(
                "no_call_record",
                Some(0),
                ".".into(),
                "G2 (rs71785313): no G2 deletion among the VCF's calls.".into(),
            );
        }
        g2(
            "not_covered",
            None,
            gq,
            "G2 (rs71785313) not read: no genotype over part of the deletion span.".into(),
        )
    }

    /// Only records overlapping the span. Like the Python window reader,
    /// stop past the window or after leaving chr22 on sorted input.
    fn window_records(
        lines: impl IntoIterator<Item = Result<String, String>>,
    ) -> Result<Vec<Vec<String>>, String> {
        let mut out = Vec::new();
        let mut seen = false;
        for line in lines {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let c = super::super::util::cells(line.strip_suffix('\n').unwrap_or(&line).as_bytes());
            if c.len() < 2 || !digits(&c[1]) {
                continue;
            }
            if c[0] != CHROM {
                if seen {
                    break;
                }
                continue;
            }
            seen = true;
            if c[1].parse::<i128>().is_ok_and(|p| p > i128::from(SPAN.1)) {
                break;
            }
            if overlaps(&c) {
                out.push(c);
            }
        }
        Ok(out)
    }

    /// The records over the G2 span, using tabix when a `.tbi` is there,
    /// otherwise streaming up to the end of the window.
    fn window(path: &str) -> Result<Vec<Vec<String>>, String> {
        let tbi = format!("{path}.tbi");
        let mut out = Vec::new();
        if std::path::Path::new(&tbi).exists() {
            let mut rd = og_io::tabix::Reader::with_index(path, &tbi).map_err(|e| e.to_string())?;
            for line in rd.query(CHROM, (SPAN.0 - 1) as u64, SPAN.1 as u64) {
                let line = line.map_err(|e| e.to_string())?;
                out.push(super::super::util::cells(&line));
            }
            return Ok(out);
        }
        let mut rd = crate::text_lines(path)?;
        window_records(std::iter::from_fn(|| {
            let mut line = String::new();
            match rd.read_line(&mut line) {
                Ok(true) => Some(Ok(line)),
                Ok(false) => None,
                Err(e) => Some(Err(format!(
                    "{}: {}",
                    base(path),
                    crate::collect::sections::gzip_message(&e)
                ))),
            }
        }))
    }

    /// G2 from the sample's gVCF when it is there, else its VCF.
    pub fn read(c: &Ctx<'_>) -> G2 {
        for (name, gvcf) in [
            (format!("vcf/{}.g.vcf.gz", c.sample), true),
            (format!("vcf/{}.vcf.gz", c.sample), false),
        ] {
            let p = c.path(&name);
            if !crate::is_file(&p) {
                continue;
            }
            return match window(&p) {
                Err(e) => G2 {
                    status: "not_read",
                    copies: None,
                    gq: ".".into(),
                    text: format!("G2 (rs71785313) not read: {name} could not be read ({e})."),
                    source: Some(name),
                    source_first: true,
                },
                Ok(recs) => G2 {
                    source: Some(name),
                    ..read_g2(&recs, gvcf)
                },
            };
        }
        G2 {
            status: "not_read",
            copies: None,
            gq: ".".into(),
            text: "G2 (rs71785313) not read: no VCF of the sample in this folder.".into(),
            source: None,
            source_first: true,
        }
    }

    /// The G1 + G2 result and its reading, from step 42's APOL1 line
    /// (`result`, `status`) and the G2 reading.
    pub fn combine(result: &str, status: &str, g2: &G2) -> (String, String) {
        let g1: Option<i64> = if status == "called" {
            match result {
                "G0/G0" => Some(0),
                "G0/G1" => Some(1),
                "G1/G1" => Some(2),
                _ => None,
            }
        } else {
            None
        };
        let two = |label: &str| format!("Two APOL1 risk copies ({label}): {TWO_COPIES}");
        let Some(g1) = g1 else {
            if g2.status == "called" && g2.copies == Some(2) {
                return (
                    "G2/G2".into(),
                    two("G2/G2")
                        + " The G1 sites were not read, but two G2 copies leave no copy for G1.",
                );
            }
            return (
                ".".into(),
                "Not read: the G1 sites of the APOL1 rule were not read.".into(),
            );
        };
        let Some(c) = g2.copies else {
            return (".".into(), "The G1 result above stands on its own.".into());
        };
        if g1 + c > 2 {
            return (
                ".".into(),
                format!(
                    "{g1} G1 and {c} G2 copies do not fit two copies of the gene (G1 and G2 never share a copy): not read. Check the calls."
                ),
            );
        }
        let mut label = std::iter::repeat_n("G0", (2 - g1 - c) as usize)
            .chain(std::iter::repeat_n("G1", g1 as usize))
            .chain(std::iter::repeat_n("G2", c as usize))
            .collect::<Vec<_>>()
            .join("/");
        let mut reading = match g1 + c {
            2 => two(&label),
            1 => format!(
                "One APOL1 risk copy ({label}): one copy did not raise the risk in Genovese 2010 (a recessive pattern)."
            ),
            _ => "No G1 and no G2 copy: no APOL1 kidney risk copy.".into(),
        };
        if g2.status == "no_call_record" {
            label += " (no G2 call)";
            reading += " No G2 deletion was called, but this VCF has no reference blocks, so a G2 copy at a position with no reads would not show.";
        }
        (label, reading)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn unindexed_window_stops_at_the_end_or_next_chromosome() {
            let call = rec(POS, REF, ALT, "0/1", ".");
            for stop in [
                format!("chr22\t{}\t.\tA\tT", SPAN.1 + 1),
                "chr23\t1\t.\tA\tT".into(),
            ] {
                let lines = [
                    Ok(call.join("\t")),
                    Ok(stop),
                    Err("must stop before reading this line".into()),
                ];
                assert_eq!(window_records(lines).unwrap(), vec![call.clone()]);
            }
        }

        #[test]
        fn unindexed_window_retains_only_overlaps() {
            use std::io::Write;
            let d = crate::testutil::Dir::new("g2-window-bounded");
            let path = d.put("calls.vcf", "##fileformat=VCFv4.2\n");
            let mut f = std::io::BufWriter::new(
                std::fs::OpenOptions::new()
                    .append(true)
                    .open(&path)
                    .unwrap(),
            );
            writeln!(f, "chr21\t1\t.\tA\tT\t.\tPASS\t.").unwrap();
            writeln!(f, "chr22\t1\t.\tA\tT\t.\tPASS\t.").unwrap();
            let block = rec(
                SPAN.0 - 20,
                "A",
                "<NON_REF>",
                "0/0",
                &format!("END={}", SPAN.1),
            );
            let call = rec(POS, REF, ALT, "0/1", ".");
            writeln!(f, "{}", block.join("\t")).unwrap();
            writeln!(f, "{}", call.join("\t")).unwrap();
            for i in 1..=250_000 {
                writeln!(f, "chr22\t{}\t.\tA\tT\t.\tPASS\t.", SPAN.1 + i).unwrap();
            }
            f.flush().unwrap();
            let records = window(&path).unwrap();
            assert_eq!(records.len(), 2);
            assert_eq!(records, vec![block, call]);
        }

        fn rec(pos: i64, r: &str, a: &str, gt: &str, info: &str) -> Vec<String> {
            [
                CHROM,
                &pos.to_string(),
                ".",
                r,
                a,
                "50",
                "PASS",
                info,
                "GT:GQ",
                &format!("{gt}:40"),
            ]
            .iter()
            .map(|x| (*x).to_owned())
            .collect()
        }

        #[test]
        fn the_marker_is_a_six_base_deletion_of_the_context() {
            let h = haplotype().expect("REF is the GRCh38 bases");
            assert_eq!(CONTEXT.len() - h.len(), 6);
        }

        #[test]
        fn g2_left_aligned_is_one_copy() {
            let g = read_g2(&[rec(POS, REF, ALT, "0/1", ".")], false);
            assert_eq!((g.status, g.copies), ("called", Some(1)));
        }

        #[test]
        fn g2_shifted_is_matched() {
            // The same deletion written three bases to the right: TAATTA>T
            // after AATAATTATAA's repeat gives the same sequence.
            let g = read_g2(&[rec(36_265_998, "ATTATA", "A", "1/1", ".")], false);
            assert_eq!(g.status, "not_read", "{g:?}");
            let shifted = rec(36_265_996, "ATAATTA", "A", "1/1", ".");
            let g = read_g2(&[shifted], false);
            assert_eq!((g.status, g.copies), ("called", Some(2)), "{g:?}");
        }

        #[test]
        fn reference_block_over_span_reads_zero() {
            let g = read_g2(
                &[rec(36_265_990, "A", "<NON_REF>", "0/0", "END=36266010")],
                true,
            );
            assert_eq!((g.status, g.copies), ("called", Some(0)));
        }

        #[test]
        fn vcf_without_record_is_no_call_record() {
            let g = read_g2(&[], false);
            assert_eq!(g.status, "no_call_record");
            let g = read_g2(&[], true);
            assert_eq!(g.status, "not_covered");
        }

        #[test]
        fn two_g2_records_are_not_read() {
            let a = rec(POS, REF, ALT, "0/1", ".");
            let b = rec(36_265_996, "ATAATTA", "A", "0/1", ".");
            let g = read_g2(&[a, b], false);
            assert_eq!(g.status, "not_read");
            assert!(g.text.contains("more than one record"));
        }

        #[test]
        fn combine_counts_both() {
            let g2 = read_g2(&[rec(POS, REF, ALT, "0/1", ".")], false);
            assert_eq!(combine("G0/G1", "called", &g2).0, "G1/G2");
            let none = read_g2(&[], false);
            let (label, reading) = combine("G0/G0", "called", &none);
            assert_eq!(label, "G0/G0 (no G2 call)");
            assert!(reading.contains("no reference blocks"));
        }
    }
}
