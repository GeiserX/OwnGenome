//! The reports show what the secondary-findings card and the carrier list
//! now say (#411: i411-display-sf, clin-33).
//!
//! The collector learned to write each mark into summary.json, but the
//! reports still printed the old rows:
//!
//! - a ClinVar hit keeps the condition ClinVar names for it, and
//!   `condition_not_acmg` when that is another disease than the one ACMG
//!   lists the gene for (a RET Hirschsprung hit in the MEN2 slot);
//! - a HIGH-impact record in TGFBR1, TGFBR2 or RBM20 carries a
//!   `mechanism_note`, and one with no GQ, DP or AD `call_unverified`;
//! - TTN and LMNA findings carry the `rule_condition` of their footnote;
//! - a recessive gene's findings carry the `pair` word, its `phase` and a
//!   `pair_note`;
//! - a run whose step 23 left no clinical VCF has its card in
//!   `sections.acmg_sf` only, built from step 6, with
//!   `high_impact_checked: false` and a `high_impact_note`;
//! - a carrier row keeps `low_confidence`, `pair` and `phase`;
//! - a two-hit line said "both copies may be hit" even when read-backed
//!   phase or a homozygous call proves it.
//!
//! A summary without the new keys keeps its old lines (the control of each
//! test), so the report-render parity goldens of such runs do not move.

use og_render::model::Summary;
use og_render::{html, text};

/// The head of a summary, then `sections`, then the tail.
fn summary(sections: &str) -> Summary {
    Summary::parse(&format!(
        r#"{{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-26 12:00 UTC",
 "run": {{"started_utc": "2026-11-26 10:00 UTC", "declared_sex": "female"}},
 "databases": {{"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null}},
 "manifest": [],
 "sections": {{ {sections} }},
 "cyp2d6": null, "not_run": [], "not_assessed": []
}}"#
    ))
    .expect("a summary")
}

/// A section of state `state` from `step`, with `data`.
fn section(key: &str, title: &str, step: &str, state: &str, data: &str) -> String {
    format!(
        r#""{key}": {{"title": "{title}", "state": "{state}", "step": "{step}", "source": "{key}/S1.tsv",
           "file_date": "2026-11-26 11:00:00", "note": null, "data": {data}}}"#
    )
}

/// Both reports, text and HTML.
fn reports(s: &Summary) -> (String, String) {
    (
        text::report(s).expect("the text report"),
        html::report(s).expect("the HTML report"),
    )
}

/// The clinical section of a run whose SF card holds `clinvar` and `high`
/// (JSON arrays of rows).
fn clinical(clinvar: &str, high: &str) -> Summary {
    summary(&section(
        "clinical",
        "Clinical Variant Filter",
        "23",
        "ok",
        &format!(
            r#"{{"variants": 3, "genes": 3, "by_impact": {{"HIGH": 1, "MODERATE": 2}},
                "acmg_sf": {{"version": "ACMG SF v3.3", "genes_on_list": 84, "rules": "assets/acmg_sf.tsv",
                  "clinvar_hits": {clinvar}, "clinvar_checked": true,
                  "high_impact": {high}, "high_impact_checked": true,
                  "high_impact_rare": 0, "high_impact_weak": 0, "high_impact_frequency_unknown": 0,
                  "structural": [], "structural_checked": false, "structural_note": null,
                  "not_reported": [], "two_hit_note": null}}}}"#
        ),
    ))
}

/// A ClinVar row of the SF card, with `extra` fields after the base ones.
fn hit(gene: &str, variant: &str, extra: &str) -> String {
    format!(
        r#"{{"gene": "{gene}", "variant": "{variant}", "genotype": "het", "significance": "Pathogenic",
            "stars": 2, "low_confidence": false, "rule": "All P and LP", "rule_note": null,
            "severity": null{extra}}}"#
    )
}

// --- the ClinVar condition and its mark ----------------------------------------------

#[test]
fn an_sf_hit_shows_the_clinvar_condition_and_the_mismatch_mark() {
    let s = clinical(
        &format!(
            "[{}]",
            hit(
                "RET",
                "chr10:43114500 C>T",
                r#", "clinvar_condition": "Hirschsprung disease", "clinvar_condition_ids": ["MedGen:C0019569"],
                   "acmg_phenotype": "familial medullary thyroid cancer, multiple endocrine neoplasia type 2",
                   "condition_not_acmg": true"#
            )
        ),
        "[]",
    );
    let (t, h) = reports(&s);
    let mark = "ClinVar Pathogenic; ClinVar condition: Hirschsprung disease, not the disease ACMG \
                lists RET for (familial medullary thyroid cancer, multiple endocrine neoplasia type \
                2): check which disease the variant was reported for";
    assert!(t.contains(mark), "{t}");
    assert!(h.contains(mark), "{h}");
}

#[test]
fn a_matching_condition_is_shown_without_the_mark() {
    let s = clinical(
        &format!(
            "[{}]",
            hit(
                "BRCA1",
                "chr17:43045712 G>A",
                r#", "clinvar_condition": "Hereditary breast ovarian cancer syndrome", "clinvar_condition_ids": [],
                   "acmg_phenotype": "hereditary breast and ovarian cancer", "condition_not_acmg": false"#
            )
        ),
        "[]",
    );
    let (t, _) = reports(&s);
    assert!(
        t.contains(
            "ClinVar Pathogenic; ClinVar condition: Hereditary breast ovarian cancer syndrome"
        ),
        "{t}"
    );
    assert!(!t.contains("not the disease ACMG lists"), "{t}");
}

/// The control: a row from before the condition keys prints as it did.
#[test]
fn an_sf_hit_without_the_new_keys_prints_as_before() {
    let s = clinical(
        &format!("[{}]", hit("BRCA1", "chr17:43045712 G>A", "")),
        "[]",
    );
    let (t, _) = reports(&s);
    assert!(
        t.contains("    BRCA1    chr17:43045712 G>A  het  ClinVar Pathogenic\n"),
        "{t}"
    );
    assert!(
        t.contains(
            "  Secondary-findings genes (ACMG SF v3.3, 84 genes): 1 ClinVar P/LP, 0 rare HIGH impact\n"
        ),
        "{t}"
    );
}

// --- the loss-of-function mark, the TTN/LMNA limit, a call with nothing to judge ------

#[test]
fn a_high_impact_row_shows_its_mechanism_note_and_an_unverified_call() {
    let high = r#"[{"gene": "TGFBR1", "variant": "chr9:99105000 C>T", "genotype": "het",
        "consequence": "stop_gained", "transcript": "MANE", "frequency": "rare", "weak_call": false,
        "call_unverified": true, "gq": null, "dp": null, "ad": null,
        "rule": "All P and LP", "rule_note": null, "severity": null,
        "mechanism_note": "TGFBR1 counts a predicted loss-of-function call with a mark: its disease is mostly caused by missense variants, so a truncating one is a less certain cause; check it with a clinician"}]"#;
    let s = summary(&section(
        "clinical",
        "Clinical Variant Filter",
        "23",
        "ok",
        &format!(
            r#"{{"variants": 1, "acmg_sf": {{"version": "ACMG SF v3.3", "genes_on_list": 84,
                "clinvar_hits": [], "high_impact": {high}, "high_impact_checked": true,
                "high_impact_rare": 0, "high_impact_weak": 0, "high_impact_frequency_unknown": 0,
                "high_impact_unverified": 1, "structural": [], "structural_checked": false,
                "not_reported": [], "two_hit_note": null}}}}"#
        ),
    ));
    let (t, h) = reports(&s);
    for page in [&t, &h] {
        assert!(
            page.contains(
                "stop_gained; no GQ, DP or AD to judge the call, not counted as rare; TGFBR1 counts \
                 a predicted loss-of-function call with a mark"
            ),
            "{page}"
        );
    }
    assert!(
        t.contains("1 with no GQ, DP or AD to judge the call)"),
        "{t}"
    );
}

#[test]
fn a_ttn_finding_shows_its_rule_condition() {
    let s = clinical(
        &format!(
            "[{}]",
            hit(
                "TTN",
                "chr2:178562000 G>A",
                r#", "rule_condition": "TTN counts only a truncating variant in its exons with high PSI; the clinical VCF holds no record of this hit, so check that it truncates""#
            )
        ),
        "[]",
    );
    let (t, h) = reports(&s);
    let shown =
        "ClinVar Pathogenic; TTN counts only a truncating variant in its exons with high PSI";
    assert!(t.contains(shown), "{t}");
    assert!(h.contains(shown), "{h}");
}

// --- the pair word of a recessive gene ---------------------------------------------------

#[test]
fn a_recessive_pair_shows_its_word_phase_and_note() {
    let pair = r#", "pair": "possibly biallelic", "phase": "phase unknown: no read-backed phase links the variants, so they may sit on the same copy (cis)",
        "pair_note": "the pair needs a HIGH-impact call of unknown frequency or with no GQ, DP or AD to judge it: check that call (its frequency in gnomAD, its reads) before you count both copies""#;
    let s = clinical(
        &format!(
            "[{}, {}]",
            hit("MUTYH", "chr1:45331556 C>T", pair),
            hit("MUTYH", "chr1:45332000 G>A", pair)
        ),
        "[]",
    );
    let (t, h) = reports(&s);
    let shown = "ClinVar Pathogenic; possibly biallelic: phase unknown: no read-backed phase links \
                 the variants, so they may sit on the same copy (cis); the pair needs a HIGH-impact \
                 call of unknown frequency";
    assert_eq!(t.matches(shown).count(), 2, "{t}");
    assert_eq!(h.matches(shown).count(), 2, "{h}");
}

// --- the card built from step 6 alone ------------------------------------------------------

/// A run whose step 23 did not run: the clinical section is not run and
/// the card lives in `sections.acmg_sf` alone.
fn step6_only(clinical_state: &str) -> Summary {
    let mut secs = vec![section(
        "acmg_sf",
        "Secondary-Findings Genes",
        "06",
        "ok",
        &format!(
            r#"{{"version": "ACMG SF v3.3", "genes_on_list": 84, "rules": "assets/acmg_sf.tsv",
                "clinvar_hits": [{}], "clinvar_checked": true, "high_impact": [],
                "high_impact_checked": false, "high_impact_rare": 0, "high_impact_weak": 0,
                "high_impact_frequency_unknown": 0, "structural": [], "structural_checked": false,
                "structural_note": null, "not_reported": [], "two_hit_note": null,
                "high_impact_note": "Rare HIGH-impact variants were not screened: step 23 (clinical filter) did not run; it needs VEP's annotation (step 13) and the VEP cache. This card holds the ClinVar P/LP hits and the class 4 and 5 structural variants only, so it can miss a predicted loss-of-function variant that is not in ClinVar, alone or as the second hit of a recessive gene."}}"#,
            hit("BRCA2", "chr13:32340300 C>T", "")
        ),
    )];
    if clinical_state != "absent" {
        secs.push(section(
            "clinical",
            "Clinical Variant Filter",
            "23",
            clinical_state,
            "{}",
        ));
    }
    summary(&secs.join(",\n"))
}

#[test]
fn a_step6_only_card_is_shown_with_its_clinvar_only_label() {
    for state in ["missing", "absent", "failed"] {
        let (t, h) = reports(&step6_only(state));
        let label = "ClinVar only, HIGH-impact not checked";
        assert!(
            t.contains(&format!(
                "  Secondary-findings genes (ACMG SF v3.3, 84 genes; {label}): 1 ClinVar P/LP\n"
            )),
            "{state}: {t}"
        );
        assert!(
            t.contains("    BRCA2    chr13:32340300 C>T  het  ClinVar Pathogenic"),
            "{state}: {t}"
        );
        assert!(
            t.contains("Rare HIGH-impact variants were not screened"),
            "{state}: {t}"
        );
        assert!(
            h.contains(&format!("Secondary-Findings Genes (ACMG SF v3.3; {label})")),
            "{state}: {h}"
        );
        assert!(
            h.contains("Rare HIGH-impact variants were not screened"),
            "{state}: {h}"
        );
    }
}

/// The control: a clinical section that holds the card shows that one, and
/// says nothing of a ClinVar-only card.
#[test]
fn the_clinical_card_wins_over_the_step6_card() {
    let (t, h) = reports(&clinical(
        &format!("[{}]", hit("BRCA1", "chr17:43045712 G>A", "")),
        "[]",
    ));
    assert!(!t.contains("ClinVar only"), "{t}");
    assert!(!h.contains("ClinVar only"), "{h}");
}

// --- the carrier list ---------------------------------------------------------------------

fn carrier(rows: &str, hom: &str) -> Summary {
    summary(&section(
        "carrier",
        "Carrier List",
        "24",
        "ok",
        &format!(
            r#"{{"genes_on_list": 3292, "list": "assets/carrier_genes.tsv", "carriers": {rows},
                "hom_or_hemizygous": {hom}, "set_aside": []}}"#
        ),
    ))
}

fn carrier_row(gene: &str, genotype: &str, extra: &str) -> String {
    format!(
        r#"{{"gene": "{gene}", "variant": "chr1:1000 C>T", "genotype": "{genotype}", "zygosity": "{genotype}",
            "significance": "Pathogenic", "stars": 2, "gq": 20, "dp": 8, "inheritance": "AR",
            "condition": "Invented condition"{extra}}}"#
    )
}

#[test]
fn a_carrier_row_shows_its_low_confidence_mark_and_its_pair() {
    let s = carrier(
        &format!(
            "[{}, {}]",
            carrier_row("GENEA", "het", r#", "low_confidence": true"#),
            carrier_row(
                "GENEB",
                "het",
                r#", "low_confidence": false, "pair": "two hits", "phase": "in trans: read-backed phase puts them on different copies (phase set 77)""#
            )
        ),
        &format!(
            "[{}]",
            carrier_row("GENEC", "hom", r#", "low_confidence": true"#)
        ),
    );
    let (t, h) = reports(&s);
    for page in [&t, &h] {
        assert!(
            page.contains(
                "GENEA (Invented condition; AR): chr1:1000 C&gt;T het, ClinVar Pathogenic; low-confidence call: check the reads"
            ) || page.contains(
                "GENEA (Invented condition; AR): chr1:1000 C>T het, ClinVar Pathogenic; low-confidence call: check the reads"
            ),
            "{page}"
        );
        assert!(
            page.contains("ClinVar Pathogenic; two hits: in trans: read-backed phase puts them on different copies (phase set 77)"),
            "{page}"
        );
        assert!(
            page.contains("two changed copies, not the pattern of a carrier; low-confidence call: check the reads"),
            "{page}"
        );
    }
}

/// The control: a carrier row with neither mark prints as it did.
#[test]
fn a_carrier_row_without_marks_prints_as_before() {
    let s = carrier(
        &format!(
            "[{}]",
            carrier_row("GENEA", "het", r#", "low_confidence": false"#)
        ),
        "[]",
    );
    let (t, _) = reports(&s);
    assert!(
        t.contains("  GENEA (Invented condition; AR): chr1:1000 C>T het, ClinVar Pathogenic\n"),
        "{t}"
    );
}

// --- the two-hit line ----------------------------------------------------------------------

#[test]
fn a_proven_two_hit_line_says_both_copies_are_hit() {
    let two = |verdict: &str| {
        format!(
            r#", "two_hits": {{"genes": [{{"gene": "GENED", "verdict": "{verdict}", "why": "a reason",
                "deletion": "a deletion", "variant": "chr1:1000 C>T", "genotype": "het", "allele": null,
                "condition": "Invented condition", "inheritance": "AR"}}],
                "deletions_read": ["read depth"], "deletions_not_read": []}}"#
        )
    };
    let with = |verdict: &str| {
        summary(&section(
            "carrier",
            "Carrier List",
            "24",
            "ok",
            &format!(
                r#"{{"genes_on_list": 3292, "carriers": [], "hom_or_hemizygous": [], "set_aside": []{}}}"#,
                two(verdict)
            ),
        ))
    };
    let (t, _) = reports(&with("biallelic"));
    assert!(
        t.contains("Biallelic: GENED (Invented condition; AR), both copies are hit: a deletion"),
        "{t}"
    );
    let (t, _) = reports(&with("possibly biallelic"));
    assert!(
        t.contains(
            "Possibly biallelic: GENED (Invented condition; AR), both copies may be hit: a deletion"
        ),
        "{t}"
    );
}
