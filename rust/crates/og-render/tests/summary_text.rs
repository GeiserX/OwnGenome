//! Report assertions from tests/test_collect_summary.py. Collect real
//! synthetic step folders, then render the same summary into both formats.

#[path = "../../og-summary/tests/port_tests_summary_support/port_tests_summary_fixtures.rs"]
mod fixtures;

use fixtures::*;
use og_render::{html, model::Summary, text};
use pycompat::json::{Dump, dumps};
use std::fs;
use std::path::Path;

fn with_reports(sample: &Sample, check: impl FnOnce(&str, &str)) {
    let value = sample.collect();
    let json = dumps(&value, &Dump::default()).expect("serialize the summary");
    let schema_path = Path::new(REPO).join("tests/schema/summary.schema.json");
    sample.put("summary.json", &json);
    let mut out = Vec::new();
    let mut err = Vec::new();
    assert_eq!(
        og_dev::schema::cli(
            &[
                schema_path.display().to_string(),
                sample.path("summary.json").display().to_string()
            ],
            &mut out,
            &mut err,
        ),
        0,
        "{}{}",
        String::from_utf8_lossy(&out),
        String::from_utf8_lossy(&err),
    );

    let summary = Summary::parse(&json).expect("parse the summary");
    let txt = text::report(&summary).expect("write the text report");
    let html = html::report(&summary).expect("write the HTML report");
    check(&txt, &html);
}

fn has(report: &str, fragments: &[&str]) {
    for fragment in fragments {
        assert!(
            report.contains(fragment),
            "missing {fragment:?} in the report"
        );
    }
}

#[test]
fn bash_and_nextflow_show_the_collected_counts_and_rows() {
    let bash = Sample::new("report-bash", "S");
    put_bash(&bash);
    with_reports(&bash, |txt, html| {
        has(
            txt,
            &[
                "Pathogenic/Likely Pathogenic hits: 3",
                "Mean depth: 29.5x",
                "Called without a function phenotype (no drug guidance): 1",
            ],
        );
        has(
            html,
            &[
                "badge-yellow\">3</span>",
                "29.5x",
                "Genes called without a function phenotype (no drug guidance)",
                "<td>PGS000018</td>",
                "<td>PGS000014</td>",
                "#chr1:20",
            ],
        );
        assert!(!html.contains("badge-green\">3</span>"));
        has(txt, &["PGS000018", "PGS000014"]);
    });

    let nf = Sample::new("report-nextflow", "S");
    put_nextflow(&nf);
    with_reports(&nf, |txt, html| {
        has(txt, &["Autosomal ROH > 5MB: 1", "31.2x"]);
        has(html, &["31.2x", "B*07:02:01", "B*08:01:01"]);
    });
}

#[test]
fn stale_and_failed_inputs_are_visible_in_both_reports() {
    let stale = Sample::new("report-stale", "S");
    put_stale(&stale);
    with_reports(&stale, |txt, html| {
        has(txt, &["[STALE: from an earlier run", "DOES NOT MATCH"]);
        has(html, &["STALE: from an earlier run"]);
    });

    let failed = Sample::new("report-failed", "S");
    put_failed(&failed);
    with_reports(&failed, |txt, html| {
        has(
            txt,
            &[
                "FAILED: step 10 failed",
                "Telomere content (TelomereHunter) (step 10 failed)",
                "failed (step 32)",
            ],
        );
        has(
            html,
            &[
                "badge-red\">Failed<",
                "step 10 failed",
                "failed (step 32)",
                "CYP2D6 Across Callers",
            ],
        );
    });
}

#[test]
fn both_reports_give_the_couple_the_same_findings_and_risk() {
    let couple = Sample::new("report-couple", "B");
    put_couple(&couple);
    with_reports(&couple, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "Partner: A",
                    "Shared: NAGA (NAGA disease): 1 in 4 per pregnancy",
                    "This sample: b-snv. A: a-snv",
                    "genetic counsellor",
                ],
            );
        }
    });
    fs::remove_file(couple.path(COUPLE)).unwrap();
    with_reports(&couple, |txt, html| {
        assert!(!txt.contains("Couple Carrier Screen"));
        assert!(!html.contains("Couple Carrier Screen"));
    });
}

#[test]
fn both_reports_show_the_read_backed_compound_het_phase() {
    for (name, nextflow, counts) in [
        (
            "report-comphet-bash",
            false,
            "1 trans, 1 cis, 1 unknown of 3 pairs",
        ),
        (
            "report-comphet-nextflow",
            true,
            "1 trans, 0 cis, 0 unknown of 1 pairs",
        ),
    ] {
        let sample = Sample::new(name, "S");
        put_comphet(&sample, nextflow);
        with_reports(&sample, |txt, html| {
            has(txt, &[counts, "G2", "19 of 20"]);
            has(html, &["G2", "19 of 20"]);
            assert!(!txt.contains("unphased single-sample data"));
        });
        if !nextflow {
            fs::remove_file(sample.path(PAIRS)).unwrap();
            with_reports(&sample, |txt, _| {
                has(txt, &["Candidates only: unphased single-sample data"]);
                assert!(!txt.contains("Phase from the reads"));
            });
        }
    }
}

#[test]
fn both_reports_explain_a_filtered_cyp2d6_call_and_the_step_36_verdict() {
    let d = Sample::new("report-cyp-unusable", "S");
    put_cyp_unusable(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "Genes called: 1/3",
                "not usable: Filter CYP2D6_depth_unreliable",
                "fewer than two usable calls",
                "Step 36: indeterminate, not passed to PharmCAT",
                "the CYP2D6 depth check failed",
            ],
        );
        has(
            html,
            &[
                "not usable: Filter CYP2D6_depth_unreliable",
                "fewer than two usable calls",
                "Not passed to PharmCAT: the CYP2D6 depth check failed",
            ],
        );
    });
    fs::remove_file(d.path(CONSENSUS)).unwrap();
    with_reports(&d, |txt, html| {
        assert!(!txt.contains("Step 36: indeterminate"));
        assert!(!html.contains("Not passed to PharmCAT: the CYP2D6 depth check failed"));
    });
}

#[test]
fn the_prs_effect_appears_only_for_an_evaluated_group() {
    let d = Sample::new("report-prs", "S");
    put_prs_evaluation(&d, "EUR", "HR 1.53 (PPM000591)");
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(report, &["HR 1.53 (PPM000591)"]);
        }
    });
    put_prs_evaluation(&d, "EAS", "");
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(report, &["evaluated in EUR"]);
            assert!(!report.contains("HR 1.53 (PPM000591)"));
        }
    });
}

#[test]
fn the_second_y_placement_and_haplocheck_have_their_own_lines() {
    let d = Sample::new("report-lineage", "S");
    put_lineage(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "Contamination (haplocheck): YES, two mtDNA haplogroups",
                "R-M269 (42 markers, QC-score 0.97)",
                "Branch: R-FT12345",
                "61 derived, 2 ancestral, 310 no call",
                "FT999",
                "FT777",
                "A0-T",
            ],
        );
        has(
            html,
            &[
                "YES, two mtDNA haplogroups",
                "Y-Chromosome Haplogroup",
                "Y Placement on the Current Public Tree",
                "<td>R-FT12345</td>",
                "<td>A0-T</td>",
            ],
        );
    });
    fs::remove_file(d.path("mito/S_haplocheck.txt")).unwrap();
    with_reports(&d, |txt, _| {
        has(txt, &["Contamination (haplocheck): not checked"])
    });
}

#[test]
fn both_reports_show_ungraded_repeats_and_a_fresh_motif_count() {
    let d = Sample::new("report-repeats", "S");
    put_ungraded_repeat(&d);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "Not graded (no thresholds in the catalog)",
                    "MUC1",
                    "SAMD12 motif count: pathogenic_motif",
                ],
            );
        }
    });
    d.old("expansion_hunter/S_str_motifs.tsv");
    with_reports(&d, |txt, html| {
        assert!(!txt.contains("SAMD12 motif count"));
        assert!(!html.contains("SAMD12 motif count"));
    });
}

#[test]
fn both_reports_render_cpsr_disagreements_and_database_versions() {
    let d = Sample::new("report-cpsr", "S");
    put_cpsr(&d);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            let report = report.replace("&#x27;", "'").replace("&gt;", ">");
            has(
                &report,
                &[
                    "BRCA2 13:32340300 G>A: final Pathogenic (ClinVar), ClinVar Pathogenic, CPSR's own Benign (BA1)",
                    "ClinVar release in CPSR's bundle: 20250309; step 6 used the ClinVar file of 20261001",
                ],
            );
        }
        has(
            txt,
            &["low-confidence call (DP 6): check the reads", "1 flagged"],
        );
    });
}

#[test]
fn export_and_mosaic_sections_validate_with_the_same_summary_schema() {
    let sample = Sample::new("report-export-mosaic", "S");
    put_export(&sample);
    put_mosaic(&sample);
    with_reports(&sample, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &["no-call in the gVCF", "a share of the cells lost the Y"],
            );
        }
    });
}

#[test]
fn both_reports_explain_construct_reads_and_always_report_calls() {
    let d = Sample::new("report-construct", "S");
    put_construct(&d);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "one copy, other copy not called ref (1/.)",
                    "known lab-construct site (APP Swedish K670N/M671L): clear",
                    "0 of 12 alt reads and 0 of 10 ref reads",
                    "chr21:25880551-25891785",
                    "ref, GQ 50 (gVCF block)",
                    "het, GQ 40, DP 30",
                    "no call: no record in the VCF or the gVCF covers the site",
                ],
            );
        }
    });
    d.put(
        CONSTRUCT_SCAN,
        "# cdna_scan: not run: no BAM\nchrom\tstart\tend\tjunctions\treads\tjunction_list\n",
    );
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(report, &["cDNA scan of every exon: not run: no BAM"]);
            assert!(!report.contains("chr21:25880551-25891785"));
        }
    });
    fs::remove_file(d.path(ALWAYS_REPORT)).unwrap();
}

#[test]
fn both_reports_put_the_smn_marker_beside_the_copy_numbers() {
    let d = Sample::new("report-parascopy", "S");
    put_parascopy(&d);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "SMN1 copies",
                    "2 (quality 33, filter PASS)",
                    "present, 5 of 40 reads carry it",
                    "SMN2 copies",
                    "does not prove 2+0",
                    "HBA1 and HBA2 together",
                ],
            );
        }
    });
    fs::remove_file(d.path(MARKER)).unwrap();
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(report, &["silent-carrier marker", "not assessed"]);
            assert!(!report.contains("present, 5 of 40 reads carry it"));
        }
    });
}

#[test]
fn the_clinical_cards_explain_which_calls_did_not_count_as_rare() {
    let d = Sample::new("report-mane", "S");
    put_clinical_mane(&d, false);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "2 rare HIGH impact",
                "Frequency unknown (not counted as rare): 1",
                "Weak genotype calls (check the reads): 1",
            ],
        );
        has(
            html,
            &[
                "Frequency unknown (not rare)",
                "Weak genotype calls",
                "GQ 12 &lt; 30",
            ],
        );
        for report in [txt, html] {
            has(report, &["frequency unknown, not counted as rare"]);
        }
    });
}

#[test]
fn the_hla_reports_distinguish_disagreement_from_no_comparison() {
    let bash = Sample::new("report-hla-qc", "S");
    put_hla_qc(&bash, false);
    with_reports(&bash, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "agree: T1K *01:01/*26:01",
                    "disagree: T1K *57:01/*08:01",
                    "not compared: T1K *07:01/*07:02",
                    "T1K gives an allele quality 0 or below",
                    "KIR2DL1*00302 (quality 60)",
                    "No allele called (absent, or too few reads): KIR2DS1",
                ],
            );
        }
    });
    let nf = Sample::new("report-hla-qc-nf", "S");
    put_hla_qc(&nf, true);
    with_reports(&nf, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "not compared: OptiType not run: no read pairs",
                    "KIR not typed: T1K found too few reads",
                ],
            );
        }
    });
}

#[test]
fn a_pypgx_only_call_is_reported_with_its_quality_and_depth_evidence() {
    let d = Sample::new("report-pgx-evidence", "S");
    put_pgx_evidence(&d, "ok");
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "WholeDel2 (pypgx, read depth)",
                    "GSTT1",
                    "GQ 12 AD 26,3",
                    "0.91 / 0.84",
                    "0.99 / 0.97",
                    "Probable (one caller: pypgx)",
                    "not passed to PharmCAT",
                ],
            );
        }
    });
    put_pgx_evidence(&d, "unreliable");
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            assert!(!report.contains("Probable (one caller: pypgx)"));
        }
    });
}

#[test]
fn the_paraphase_report_distinguishes_disagreement_agreement_and_one_caller() {
    let d = Sample::new("report-paraphase", "S");
    put_paraphase(&d, 1);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "SMN1 by both tools",
                    "they DISAGREE",
                    "HBA2 lost",
                    "Variants on phased haplotypes",
                ],
            );
        }
    });
    put_paraphase(&d, 2);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(report, &["the two agree"]);
        }
    });
    fs::remove_file(d.path(SMN)).unwrap();
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            assert!(!report.contains("SMN1 by both tools"));
        }
    });
}

#[test]
fn repeat_reports_show_the_sizing_limits_and_extra_tract_readings() {
    let d = Sample::new("report-repeat-sizing", "S");
    put_repeat_sizing(&d);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "size not determined (flanking reads only)",
                    "In an intermediate range",
                    "HTT, DMPK",
                    "RFC1 motif count: benign_motif",
                ],
            );
        }
    });

    let d = Sample::new("report-extra-repeats", "S");
    put_extra_repeats(&d, true);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &["TG11/TG12 5T/7T", "5T: 1 copy", "*1/*28", "5-HTTLPR"],
            );
        }
    });
    put_extra_repeats(&d, false);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(report, &["TG11-7T/TG11-7T"]);
            assert!(!report.contains("5T: 1 copy"));
        }
    });
}
