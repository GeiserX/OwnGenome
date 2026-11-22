//! Report assertions from tests/test_collect_summary.py. Collect real
//! synthetic step folders, then render the same summary into both formats.

#[path = "../../og-summary/tests/port_tests_summary_support/port_tests_summary_fixtures.rs"]
mod fixtures;

use fixtures::*;
use og_render::{html, model::Summary, text};
use pycompat::json::{Dump, Value, dumps, loads};
use std::fs;
use std::path::Path;

fn with_reports(sample: &Sample, check: impl FnOnce(&str, &str)) {
    render_value(sample, &sample.collect(), check);
}

fn render_value(sample: &Sample, value: &Value, check: impl FnOnce(&str, &str)) {
    let json = dumps(value, &Dump::default()).expect("serialize the summary");
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

    render_reports(&json, check);
}

fn render_reports(json: &str, check: impl FnOnce(&str, &str)) {
    let summary = Summary::parse(json).expect("parse the summary");
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
                "Genes called without a function phenotype (no drug guidance)</span><span class=\"value\">1<",
                "<td>PGS000018</td>",
                "<td>PGS000014</td>",
                "#chr1:20",
            ],
        );
        assert!(!html.contains("badge-green\">3</span>"));
        let rows: Vec<_> = html
            .lines()
            .filter(|line| line.trim_start().starts_with("<tr><td>chr"))
            .collect();
        assert_eq!(rows.len(), 3);
        assert!(rows[0].trim_end().ends_with("<td>3</td><td></td><td>gnomAD AF unknown</td><td></td><td><a class=\"locus\" href=\"http://127.0.0.1:8765/#chr1:20\" target=\"_blank\" rel=\"noopener noreferrer\">chr1:20</a></td></tr>"));
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
        has(html, &["<div class=\"stale\">STALE: from an earlier run"]);
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
        let telomere = card(html, "Telomere content (relative)");
        has(telomere, &["badge-red\">Failed<", "step 10 failed"]);
        assert!(!telomere.contains("Not run"));
        let callers = card(html, "CYP2D6 Across Callers");
        has(callers, &["failed (step 32)", "FAILED: step 32 failed"]);
        // The CYP2D6 section places its status badge on its first card.
        has(
            card(html, "pypgx Pharmacogenomics"),
            &["badge-red\">Failed<"],
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
            has(
                txt,
                &[
                    &format!("  Phase from the reads: {counts}\n"),
                    "G2",
                    "19 of 20",
                ],
            );
            has(
                html,
                &["<td>G2</td><td>chr1:50 A&gt;G + chr1:90 C&gt;T</td><td>19 of 20</td>"],
            );
            if !nextflow {
                has(
                    txt,
                    &[
                        "G2       chr1:50 A>G + chr1:90 C>T  in trans on 19 of 20 reads or read pairs",
                    ],
                );
                has(html, &["1 / 1 / 1"]);
            }
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
    put_prs_evaluation(&d, "EUR", PRS_EUR_EFFECT);
    with_reports(&d, |txt, html| {
        has(txt, &[&format!("effect {PRS_EUR_EFFECT}")]);
        has(
            html,
            &[
                "<th>Effect in your group</th>",
                &format!("<td>{PRS_EUR_EFFECT}</td>"),
            ],
        );
    });
    put_prs_evaluation(&d, "EAS", "");
    let missing = "no published OR or HR for EAS; evaluated in EUR, AFR, AMR";
    with_reports(&d, |txt, html| {
        has(txt, &[missing]);
        has(html, &[&format!("<td>{missing}</td>")]);
        for report in [txt, html] {
            assert!(!report.contains("HR 1.53"));
        }
    });
    put_prs_without_evaluations(&d);
    with_reports(&d, |txt, html| {
        has(txt, &["percentile 81.2 (EUR)"]);
        assert!(!txt.contains("Effect: the odds ratio"));
        assert!(!html.contains("Effect in your group"));
        for report in [txt, html] {
            assert!(!report.contains("no published OR"));
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
                "Branch: R-FT12345 (level 24, YFull YTree 14.06.0)",
                "61 derived, 2 ancestral, 310 no call",
                "Ancestral inside the terminal branch: FT999",
                "outvoted by the SNPs below it: FT777",
                "A0-T",
            ],
        );
        has(
            html,
            &[
                "YES, two mtDNA haplogroups",
                "Y Placement on the Current Public Tree",
                "<td>R-FT12345</td>",
                "<td>A0-T</td>",
            ],
        );
        has(
            card(html, "Y-Chromosome Haplogroup"),
            &["<span class=\"label\">Haplogroup</span><span class=\"value\">R-M269</span>"],
        );
    });
    put_unplaced_y(&d);
    with_reports(&d, |txt, _| {
        has(
            txt,
            &[
                "Haplogroup: insufficient markers (3 markers",
                "no branch has a derived call",
            ],
        );
        assert!(!txt.contains("Ancestral inside the terminal branch"));
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
        let why = "the reference repeat is 1200 bases, longer than 300: short reads cannot size it";
        has(report_row(txt, "MUC1"), &["not graded", why]);
        has(
            txt,
            &[
                &format!("Not graded (no thresholds in the catalog): 1 (MUC1: {why})"),
                "SAMD12 motif count: pathogenic_motif: motif TGAAA (pathogenic) in 5 of 5 repeat reads",
                og_render::notes::NOT_GRADED_NOTE,
            ],
        );
        has(
            card(html, "Repeat Expansions"),
            &[
                "Not graded (no thresholds in the catalog)",
                "<td>MUC1</td><td>20/25</td><td>not graded</td>",
                "SAMD12 motif count: pathogenic_motif",
            ],
        );
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
            &[
                "Rows where ClinVar and CPSR's own class disagree on P/LP: 3",
                "low-confidence call (DP 6): check the reads",
                "cDNA read check of the P/LP rows (alt reads clipped at exon edges): 1 flagged",
            ],
        );
    });
    d.gzip("cpsr/S.cpsr.grch38.classification.tsv.gz", "VAR_ID\tSYMBOL\tCLASSIFICATION\n13_32340300_G_A\tBRCA2\tPathogenic\n5_112800000_T_C\tAPC\tBenign\n");
    with_reports(&d, |txt, _| {
        has(txt, &["has no CPSR_CLASSIFICATION column"])
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
        has(
            txt,
            &[
                "Known deep intronic sites (always reported)",
                "one copy, other copy not called ref (1/.)  Pathogenic",
                "reads: known lab-construct site (APP Swedish K670N/M671L): clear",
                "cDNA scan of every exon: 1 loci with reads across two or more exon junctions (read 3 exon regions on 1 contigs; 1 loci listed)",
                "chr21:25880551-25891785: 3 junctions, 24 reads",
                og_render::notes::CDNA_SCAN_NOTE,
            ],
        );
        assert!(!txt.contains("chr7:117548628 GTTT>G  het"));
        assert_eq!(html.matches("<td>NM_").count(), 3);
        has(
            html,
            &[
                "<td>one copy, other copy not called ref (1/.)</td>",
                "<td>chr21:25880551-25891785</td><td>3</td><td>24</td>",
                "Hits at a known lab-construct site (check the reads)</span><span class=\"value\"><span class=\"badge badge-yellow\">1<",
            ],
        );
        for (gene, call) in [
            ("USH2A", "ref, GQ 50 (gVCF block)"),
            ("CFTR", "het, GQ 40, DP 30"),
            (
                "GLA",
                "no call: no record in the VCF or the gVCF covers the site",
            ),
        ] {
            assert!(
                txt.lines()
                    .any(|line| line.contains(gene) && line.contains(call))
            );
            let row = format!("<td>{gene}</td>");
            let call_cell = format!("<td>{call}</td>");
            assert!(
                html.split("<tr>").any(|cells| {
                    let cells = cells.split("</tr>").next().unwrap();
                    cells.contains(&row) && cells.contains(&call_cell)
                }),
                "missing {call:?} in the {gene} row"
            );
        }
        for report in [txt, html] {
            has(
                report,
                &[
                    "one copy, other copy not called ref (1/.)",
                    "known lab-construct site (APP Swedish K670N/M671L): clear",
                    "0 of 12 alt reads and 0 of 10 ref reads",
                    "chr21:25880551-25891785",
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
            assert!(!report.contains("Loci with spliced reads"));
        }
    });
    fs::remove_file(d.path(ALWAYS_REPORT)).unwrap();
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            assert!(!report.contains("Known deep intronic sites (always reported)"));
        }
    });
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
                "0 ClinVar P/LP, 2 rare HIGH impact (listed apart: 2 weak call(s), 1 of unknown frequency)",
                "Frequency unknown (not counted as rare): 1",
                "Weak genotype calls (check the reads): 1",
                "Frequency unknown (prioritized, not counted as rare): 1",
                "weak call (GQ 12 < 30): check the reads",
            ],
        );
        has(
            html,
            &[
                "Frequency unknown (not rare)",
                "Weak genotype calls",
                "weak call (GQ 12 &lt; 30): check the reads",
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
        has(
            txt,
            &[
                "OptiType check (HLA_QC=true, two fields):",
                "HLA-A    agree: T1K *01:01/*26:01, OptiType *26:01/*01:01",
                "HLA-B    disagree: T1K *57:01/*08:01, OptiType *58:01/*08:01 (the callers disagree at two fields; step 36 does not pass it to PharmCAT)",
                "HLA-C    not compared: T1K *07:01/*07:02, OptiType *07:01/*07:02 (T1K gives an allele quality 0 or below)",
                "KIR genes (T1K, KIR=true):",
                "KIR2DL1  KIR2DL1*00302 (quality 60) / KIR2DL1*0010101 (quality 35)",
                "KIR3DL1  KIR3DL1*0150201 (quality 50)",
                "No allele called (absent, or too few reads): KIR2DS1",
                og_render::notes::KIR_NOTE,
            ],
        );
        has(
            card(html, "HLA Typing"),
            &[
                "HLA-A: T1K vs OptiType",
                "agree: T1K *01:01/*26:01, OptiType *26:01/*01:01",
                "HLA-B: T1K vs OptiType",
                "disagree: T1K *57:01/*08:01, OptiType *58:01/*08:01",
                "not compared: T1K *07:01/*07:02",
                "T1K gives an allele quality 0 or below",
                "KIR2DL1*00302 (quality 60) / KIR2DL1*0010101 (quality 35)",
                "No allele called (absent, or too few reads): KIR2DS1",
                og_render::notes::KIR_NOTE,
            ],
        );
    });

    let nf = Sample::new("report-hla-qc-nf", "S");
    put_hla_qc(&nf, true);
    with_reports(&nf, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "not compared: OptiType not run: no read pairs at HLA-A, HLA-B or HLA-C in this BAM",
                    "KIR not typed: T1K found too few reads at the KIR genes",
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
        let block = txt
            .split_once("## PGx copy number")
            .unwrap()
            .1
            .split("##")
            .next()
            .unwrap();
        has(block, &["GSTM1    WholeDel2 (pypgx, read depth)"]);
        assert!(!block.contains("UGT2B17"));
        has(html, &["PGx Copy Number", "WholeDel2 (pypgx, read depth)"]);
        assert!(!html.contains("CNVpytor deletion chr4:68500001-68650000"));
        let depth = "CYP2D6 depth check: ok;";
        let shares = "share of the depth at MAPQ >= 1 / >= 20: Cyrius's region with REP6 0.91 / 0.84, gene body without REP6 0.99 / 0.97";
        let alt = "Not assessable on this reference: GSTT1 (GRCh38 has it only on the ALT contig chr22_KI270879v1_alt";
        let weak = "Low-confidence: CYP2C19 *1/*2 rests on chr10:94781859 G>A GQ 12 AD 26,3: check the reads";
        let text_depth = txt
            .lines()
            .find(|line| line.trim_start().starts_with(depth))
            .unwrap();
        has(text_depth, &[shares]);
        let html_depth = html.lines().find(|line| line.contains(depth)).unwrap();
        has(html_depth, &[&og_render::escape::e(shares)]);
        let probable = "*1/*2, Normal Metabolizer; not passed to PharmCAT";
        has(
            txt,
            &[
                depth,
                alt,
                weak,
                &format!("Probable (one caller: pypgx): {probable}"),
            ],
        );
        has(
            html,
            &[
                &og_render::escape::e(depth),
                &og_render::escape::e(alt),
                &og_render::escape::e(weak),
                "Probable (one caller: pypgx)",
                probable,
            ],
        );
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
        has(
            txt,
            &[
                "\n  HBA2: 1 (called); genotype -a/aa (3 alpha-globin genes); 4.2 kb deletion (-a4.2): HBA2 lost\n",
                "\n  SMN1 by both tools: Paraphase 2, Parascopy 1 (quality 33, filter PASS): they DISAGREE: check both before using either\n",
                "\n  Variants on phased haplotypes: 1\n",
            ],
        );
        has(
            card(html, "Paralog Genes (Paraphase)"),
            &[
                "<td>HBA2</td><td>1</td>",
                "they DISAGREE",
                "PacBio HiFi reads only",
                "HBA2 lost",
                "Variants on phased haplotypes",
            ],
        );
    });
    put_paraphase(&d, 2);
    with_reports(&d, |txt, html| {
        let reading = "Parascopy 2 (quality 33, filter PASS): the two agree";
        has(
            txt,
            &[&format!("\n  SMN1 by both tools: Paraphase 2, {reading}\n")],
        );
        has(card(html, "Paralog Genes (Paraphase)"), &[reading]);
    });
    fs::remove_file(d.path(SMN)).unwrap();
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            assert!(!report.contains("SMN1 by both tools"));
        }
    });
    let missing = Sample::new("report-paraphase-missing", "S");
    with_reports(&missing, |txt, _| {
        has(txt, &["\n  - Paralog genes (Paraphase, PacBio HiFi)\n"]);
    });
}

fn card<'a>(html: &'a str, title: &str) -> &'a str {
    html.split_once(&format!("<h2>{title}</h2>"))
        .unwrap_or_else(|| panic!("missing {title} card"))
        .1
        .split_once("\n  </div>")
        .unwrap_or_else(|| panic!("unclosed {title} card"))
        .0
}

fn report_row<'a>(text: &'a str, locus: &str) -> &'a str {
    let mut rows = text
        .lines()
        .filter(|line| line.trim_start().starts_with(&format!("{locus} ")));
    let row = rows.next().unwrap_or_else(|| panic!("missing {locus} row"));
    assert!(rows.next().is_none(), "more than one {locus} row");
    row
}

#[test]
fn stranger_reports_distinguish_fresh_outdated_and_stale_results() {
    let d = Sample::new("report-stranger-status", "S");
    put_stranger(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &["LOCUS_B  20/75  full_mutation", "can be wrong at some loci"],
        );
        has(
            card(html, "Repeat Expansions"),
            &[
                "<td>LOCUS_B</td><td>20/75</td><td>full_mutation</td>",
                "can be wrong at some loci",
            ],
        );
    });

    d.date(STRANGER, STARTED - 120);
    with_reports(&d, |txt, html| {
        let html = card(html, "Repeat Expansions");
        has(txt, &["rerun step 9b"]);
        has(html, &["rerun step 9b"]);
        assert!(!txt.contains("did not run"));
        assert!(!html.contains("not run"));
    });

    d.put(
        RUN,
        format!("meta\tstarted_epoch\t{STARTED}\nstep\t09\tskipped\n"),
    );
    with_reports(&d, |_, html| {
        let html = card(html, "Repeat Expansions");
        has(html, &["Stale"]);
        assert!(!html.contains("Complete"));
    });
}

#[test]
fn repeat_reports_show_the_sizing_limits_and_extra_tract_readings() {
    let d = Sample::new("report-repeat-sizing", "S");
    put_repeat_sizing(&d);
    with_reports(&d, |txt, html| {
        let html = card(html, "Repeat Expansions");
        has(
            txt,
            &[
                "In an intermediate range (bin/str_intermediate_ranges.tsv): HTT, DMPK",
                "RFC1 motif count: benign_motif: motif AAAAG (benign) in 22 of 24 repeat reads (92%): the flag is not supported",
                ">=80/>=130",
            ],
        );
        has(
            html,
            &[
                "size not determined (flanking reads only)",
                "In an intermediate range",
                "HTT, DMPK",
                "RFC1 motif count: benign_motif",
                "<td>&gt;=80/&gt;=130</td>",
            ],
        );
        assert_eq!(html.matches("<tr><td>").count(), 10);
        let flanking = report_row(txt, "LOCUS_FL");
        assert_ne!(flanking.split_whitespace().nth(2), Some("normal"));
        has(
            flanking,
            &["not sized  size not determined (flanking reads only)"],
        );
        let spanning = report_row(txt, "LOCUS_SP");
        assert_eq!(spanning.split_whitespace().nth(2), Some("normal"));
        assert!(!spanning.contains("size not determined"));
        has(
            report_row(txt, "HTT"),
            &[
                "27-35: intermediate (not disease-causing; may expand in offspring)",
                "https://www.ncbi.nlm.nih.gov/books/NBK1305/",
            ],
        );
    });

    let d = Sample::new("report-extra-repeats-exa", "S");
    put_extra_repeats(&d, extra_repeats_exa());
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "CFTR intron 9 (legacy IVS8) TG and poly-T tract: TG11/TG12 5T/7T",
                "    5T: 1 copy",
                "UGT1A1 promoter TA repeat: *1/*28",
                "SLC6A4 5-HTTLPR promoter length polymorphism: not tested",
            ],
        );
        has(
            card(html, "Repeat Expansions"),
            &["TG11/TG12 5T/7T", "5T: 1 copy", "*1/*28", "5-HTTLPR"],
        );
    });

    let d = Sample::new("report-extra-repeats-exb", "S");
    put_extra_repeats(&d, extra_repeats_exb());
    with_reports(&d, |txt, _| assert!(!txt.contains("5T:")));
}

#[test]
fn empty_repeat_vcf_lists_unsized_extra_tracts_and_missing_key_loci() {
    let d = Sample::new("report-empty-repeats", "S");
    put_extra_repeats(&d, ExtraRepeats::default());
    with_reports(&d, |txt, _| {
        for locus in EH_LOCI {
            has(txt, &[&format!("{locus:<8} not in output")]);
        }
        has(
            txt,
            &["CFTR intron 9 (legacy IVS8) TG and poly-T tract: not sized"],
        );
    });
}

#[test]
fn clinvar_reports_show_call_quality_and_allele_frequency() {
    let d = Sample::new("report-clinvar-quality", "S");
    put_clinvar_quality(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "frequency: gnomAD AF 0.03; common: unlikely to be a high-penetrance finding",
                "frequency: gnomAD AF unknown",
                "call: low-confidence call (GQ 12, DP 6): check the reads",
                "Low-confidence calls (GQ below 30 or DP below 10), listed but not counted above: 1",
            ],
        );
        has(
            html,
            &[
                "<td>gnomAD AF 0.03; common: unlikely to be a high-penetrance finding</td>",
                "<th>gnomAD frequency</th>",
                "<td>low-confidence call (GQ 12, DP 6): check the reads</td>",
            ],
        );
        for report in [txt, html] {
            has(
                report,
                &["a heterozygous hit in a dominant gene as well as a homozygous one"],
            );
        }
    });
    fs::remove_dir_all(d.path("vep")).unwrap();
    with_reports(&d, |txt, _| {
        has(
            txt,
            &["gnomAD AF: unknown for every hit (no VEP output (step 13) in this folder)"],
        )
    });
}

#[test]
fn mito_reports_show_polyc_and_annotation_status() {
    let d = Sample::new("report-mito-annotation", "S");
    put_mito_calls(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "Heteroplasmic (AF 0.05 to 0.95, outside the D-loop poly-C tracts): 1",
                "Length heteroplasmy in the D-loop poly-C tracts (m.303-315, m.16184-16193; common, usually technical): 2",
                "ClinVar and gnomAD: not annotated (no annotation table",
                "MITO_VARIANTS",
            ],
        );
        has(html, &["Length heteroplasmy in the D-loop poly-C tracts"]);
    });
    put_mito_databases(&d);
    annotate_mito(&d, true);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "ClinVar pathogenic or likely pathogenic: 1",
                "chrM:3243 A>G  AF 0.3  ClinVar Pathogenic (reviewed_by_expert_panel; MELAS)  gnomAD AF hom 0.00002, het 0.0003",
                "2026-10-05",
            ],
        );
        assert!(!txt.contains("fileDate 2026-10-05"));
        has(html, &["<td>MELAS</td>"]);
    });
    annotate_mito(&d, false);
    with_reports(&d, |txt, _| {
        has(
            txt,
            &["ClinVar not annotated: no ", "gnomAD not annotated: no "],
        )
    });
    d.date(MITO_ANNOTATION, STARTED - 120);
    with_reports(&d, |txt, _| {
        has(txt, &["older than the chrM calls: rerun step 20"])
    });
}

#[test]
fn parascopy_catalog_reports_call_readings_and_unresolved_regions() {
    let d = Sample::new("report-parascopy-catalog", "S");
    put_parascopy_catalog(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "HBA1 and HBA2 together (alpha-globin genes, from depth): 4",
                "  GBA1, GBAP1 (agCN): 4 (2 region(s), quality 38 or more, filter PASS). Reading: agCN of GBA1",
                "  AMY1A, AMY1B, AMY1C (agCN): not assessed: no read over",
            ],
        );
        let c4 = report_row(txt, "C4A,");
        has(c4, &["C4A, C4B (agCN): chr6"]);
        assert!(!c4.contains("Reading: "));
        has(
            card(html, "SMN1/SMN2 Copy Number (Parascopy)"),
            &[
                "<td>GBA1, GBAP1</td>",
                "<th>Copies together (agCN)</th>",
                "Gaucher",
            ],
        );
    });
}

#[test]
fn paralog_reports_explain_missing_failed_and_paraphase_only_steps() {
    let d = Sample::new("report-paralog-availability", "S");
    d.gzip(
        "vcf/S.vcf.gz",
        &format!("{VCF}chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"),
    );
    with_reports(&d, |txt, html| {
        let line = "Paralog copy number: not checked (neither step 35 nor step 35b ran)";
        has(txt, &[&format!("## Paralog Copy Number\n---\n  {line}\n")]);
        has(card(html, "Paralog Copy Number"), &["Not checked", line]);
    });
    d.put(
        RUN,
        format!("meta\tstarted_epoch\t{STARTED}\nstep\t35\tfailed (exit 1)\n"),
    );
    with_reports(&d, |txt, _| {
        has(
            txt,
            &["Paralog copy number: not checked (step 35 failed, step 35b not run)"],
        )
    });
    fs::remove_file(d.path(RUN)).unwrap();

    d.put(PARAPHASE_CN, format!("{PARAPHASE_HEADER}{HBA1}{HBA2}"));
    with_reports(&d, |txt, html| {
        has(txt, &["HBA2: 2 (called)"]);
        assert!(!txt.contains("not checked"));
        assert!(!html.contains("Not checked"));
    });
    d.put(
        PARAPHASE_CN,
        format!("{PARAPHASE_HEADER}smn1\tSMN1\t2\tcalled\t.\t4\t90\t31\n"),
    );
    with_reports(&d, |txt, _| {
        assert!(!txt.contains("Paralog copy number: not checked"))
    });
}

#[test]
fn paralog_reports_explain_hba_marker_and_copy_number_readings() {
    let d = Sample::new("report-paralog-readings", "S");
    d.put(SMN, format!("{SMN_HEADER}{SMN_ELSEWHERE}{SMN_EXON_7}"));
    let after_cn = "  SMN1 copies: 2 (quality 33, filter PASS)\n  SMN1 g.27134T>G silent-carrier marker (rs143838139): ";
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "## SMN1/SMN2 Copy Number (Parascopy, short reads)",
                "  SMN2 copies: 2 (quality 30, filter PASS)\n",
                "  SMN1 and SMN2 together (agCN) at exon 7: 4 (quality 60, filter PASS)\n",
                "  Reading: two SMN1 copies: not a carrier",
                &format!("{after_cn}not assessed\n"),
                "step 35 wrote no S_smn1_marker.tsv",
                "SMN1 marker caveat: The marker raises the chance that two SMN1 copies sit on one",
                "HBA1/HBA2 copy number (alpha-thalassaemia deletions): neither step 35's read-depth count",
            ],
        );
        assert!(!txt.contains("HBA1 and HBA2 together"));
        assert!(!txt.contains("not checked"));
        assert!(!html.contains("Not checked"));
        has(
            card(html, "SMN1/SMN2 Copy Number (Parascopy)"),
            &[
                ">SMN1 copies<",
                "2 (quality 33, filter PASS)",
                ">SMN2 copies<",
                "not a carrier by copy number",
                "CYP21A2 and the gene-dosage loci have no Parascopy model",
                "step 35 counts them from read depth against flanking sequence",
            ],
        );
    });
    put_parascopy(&d);
    with_reports(&d, |txt, html| {
        let caveat = og_render::sections::smn1::SMN1_MARKER_CAVEAT;
        has(
            txt,
            &[
                &format!("{after_cn}present, 5 of 40 reads carry it, about 0.5 copies\n"),
                "  SMN1 marker reading: SMN1 g.27134T>G (rs143838139): present",
                &format!("\n  SMN1 marker caveat: {caveat}\n"),
                "  HBA1 and HBA2 together (alpha-globin genes, from depth): 4 (four alpha-globin genes (aa/aa): no deletion seen; estimate 4.02)\n",
            ],
        );
        has(
            card(html, "SMN1/SMN2 Copy Number (Parascopy)"),
            &[
                ">SMN1 g.27134T&gt;G silent-carrier marker (rs143838139)<",
                "present, 5 of 40 reads carry it",
                &og_render::escape::e(caveat),
                "does not prove 2+0",
                ">HBA1 and HBA2 together (alpha-globin genes, from depth)<",
                "estimate 4.02",
            ],
        );
    });
    d.put(MARKER, "metric\tvalue\nstatus\tmaybe\n");
    with_reports(&d, |txt, _| {
        has(
            txt,
            &[
                &format!("{after_cn}not assessed\n"),
                "S_smn1_marker.tsv has no status row",
            ],
        )
    });
    d.put(HBA, HBA_DEPTH_3);
    with_reports(&d, |txt, _| {
        has(
            txt,
            &[
                "from depth): 3 (three alpha-globin genes: one gene lost on one chromosome",
                "a silent carrier; estimate 3.01)",
            ],
        )
    });
    d.put(HBA, HBA_LOW);
    with_reports(&d, |txt, _| {
        has(txt, &["from depth): not assessed: flank depth under 8x"])
    });
    d.put(SMN, format!("{SMN_HEADER}{SMN_CARRIER}"));
    with_reports(&d, |txt, _| {
        has(
            txt,
            &["  Reading: one SMN1 copy: the pattern of an SMA carrier\n"],
        )
    });
}

#[test]
fn copy_number_reports_distinguish_dosage_interval_and_missing_tables() {
    let only = Sample::new("report-dosage-only", "S");
    put_dosage_only(&only);
    with_reports(&only, |txt, html| {
        has(
            txt,
            &["HP (haptoglobin) copies (from depth): 3 (estimate 3.04)"],
        );
        has(html, &[">HP (haptoglobin) copies (from depth)<"]);
    });

    let d = Sample::new("report-cn-tables", "S");
    d.put(SMN, format!("{SMN_HEADER}{SMN_EXON_7}"));
    with_reports(&d, |txt, _| {
        has(
            txt,
            &[
                "  CYP21A2 (RCCX module, from depth and six PSV sites): not assessed (step 35's <id>_cyp21a2_psv.tsv is not in this run's output)\n",
                "  Gene-dosage loci (HP, LCE3B/LCE3C, CCL3L1): not assessed (step 35's <id>_cn_dosage.tsv",
                "  Copy-number intervals (ClinGen recurrent regions, carrier deletions, AZF; from depth): not assessed (step 35's <id>_cn_intervals.tsv",
            ],
        );
        assert!(!txt.contains("is not assessed"));
        assert!(!txt.contains("CYP21A2 has no Parascopy model"));
    });
    put_cn_tables(&d);
    with_reports(&d, |txt, html| {
        let html = card(html, "SMN1/SMN2 Copy Number (Parascopy)");
        has(
            txt,
            &[
                "  CYP21A2 (RCCX module, from depth and six PSV sites): 1 CYP21A2 copies, module total 3 of 4, finding deletion. Reading: one RCCX module fewer",
                "  HP (haptoglobin) copies (from depth): 3 (estimate 3.04). Reading: three copies of the duplicated HP exons: Hp2-1\n",
                "  LCE3B/LCE3C copies (from depth): 1 (estimate 0.98). Reading: LCE3B and LCE3C deleted",
                "  CCL3L1/CCL3L3 copies (from depth): not assessed (check_failed)\n",
                "  Copy-number intervals (ClinGen recurrent regions, carrier deletions, AZF; from depth): 3 rows (2 called, 1 not assessed); 1 with a deletion or a duplication\n",
                "    Deletion: 22q11.2 recurrent (DGS/VCFS) region (proximal, A-D) (includes TBX1) (22q11.2_LCR_A-LCR_D, chr22:18924718-21111383), 1 of 2 copies\n",
            ],
        );
        assert!(!txt.contains("CFTRdele2,3"));
        has(
            html,
            &[
                ">CYP21A2 (RCCX module, from depth and six PSV sites)<",
                "1 CYP21A2 copies, module total 3 of 4, finding deletion",
                ">HP (haptoglobin) copies (from depth)<",
                ">CCL3L1/CCL3L3 copies (from depth)<",
                "3 rows (2 called, 1 not assessed); 1 with a deletion or a duplication",
                "<td>22q11.2_LCR_A-LCR_D</td><td>chr22:18924718-21111383</td><td>1 of 2</td><td>deletion</td>",
            ],
        );
    });
    d.put(
        CYP21A2,
        "metric\tvalue\nstatus\tnot assessed\ncyp21a2_copies\t.\nreason\tpsv_low_depth\n",
    );
    d.put(INTERVALS, format!("{IV_COLS}{INTERVAL_UNCHANGED}"));
    with_reports(&d, |txt, _| {
        has(
            txt,
            &[
                "  CYP21A2 (RCCX module, from depth and six PSV sites): not assessed (psv_low_depth)\n",
                "from depth): 1 rows (1 called); no deletion or duplication among the 1 called rows\n",
            ],
        )
    });
    d.put(INTERVALS, format!("{IV_COLS}CFTR_dele2_3\tCFTRdele2,3\tchr7\t117479000\t117481000\tnot assessed\t.\t2\t.\t.\tcheck_failed\n\
        AZFc\tAZFc\tchrY\t24000000\t26000000\tnot assessed\t.\t1\t.\t.\tcheck_failed\n"));
    with_reports(&d, |txt, html| {
        has(
            txt,
            &["from depth): 2 rows (2 not assessed); none called\n"],
        );
        has(
            card(html, "SMN1/SMN2 Copy Number (Parascopy)"),
            &["2 rows (2 not assessed); none called"],
        );
        assert!(!txt.contains("no deletion or duplication"));
    });
}

#[test]
fn combined_alleles_are_reported_beside_their_clinvar_stars() {
    let d = Sample::new("report-combined-alleles", "S");
    d.put(CLINVAR, combined_hits("0/1"));
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                &format!(
                    "(part of a combined allele: chr22:30000038 GCCAT>CATGGG{COMBINED_UNPHASED})"
                ),
                &format!(
                    "(combined allele of chr22:30000037 AGC>A + chr22:30000042 T>TGGG{COMBINED_UNPHASED})"
                ),
                "possibly one allele with chr22:30000305 C>T: check ClinVar for a combined record",
            ],
        );
        has(
            html,
            &[
                "<th>Stars</th><th>Combined allele</th>",
                &format!(
                    "<td>part of a combined allele: chr22:30000038 GCCAT&gt;CATGGG{COMBINED_UNPHASED}</td>"
                ),
            ],
        );
        let snv = html
            .lines()
            .find(|line| {
                line.trim_start().starts_with("<tr><td>chr22") && line.contains("<td>30000210</td>")
            })
            .unwrap();
        assert!(snv.trim_end().ends_with("<td></td><td><a class=\"locus\" href=\"http://127.0.0.1:8765/#chr22:30000210\" target=\"_blank\" rel=\"noopener noreferrer\">chr22:30000210</a></td></tr>"));
    });
}

#[test]
fn hla_reports_show_each_alleles_quality_without_unrun_opt_ins() {
    let d = Sample::new("report-hla-allele-quality", "S");
    put_hla_allele_quality(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "A*01:01:01 (quality 60)\n",
                "B*57:01:01 (quality 0) / B*08:01:01 (quality 40); low confidence",
                "C*07:02:01 (quality 30); low confidence (T1K quality 0 or below)\n",
            ],
        );
        assert_eq!(
            txt.matches("low confidence (T1K quality 0 or below), not passed to PharmCAT")
                .count(),
            1
        );
        has(
            html,
            &[
                "B*57:01:01 (quality 0)",
                "HLA typing from short-read WGS is approximate",
            ],
        );
        assert_eq!(html.matches("not passed to PharmCAT</span>").count(), 1);
        for report in [txt, html] {
            assert!(!report.contains("OptiType"));
            assert!(!report.contains("KIR"));
        }
    });
}

#[test]
fn roh_and_telomere_reports_explain_masking_and_the_read_floor() {
    let d = Sample::new("report-roh-and-telomere", "S");
    put_masked_roh(&d);
    put_telomere(&d, 40);
    let warning = "too few telomeric reads for a stable estimate";
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "Left out in centromere and heterochromatin bands: 46.5 MB",
                "No population allele frequencies were used",
                "Segments on chrX/chrY, not counted: 1",
                "Intratelomeric reads: 40",
                warning,
            ],
        );
        has(html, &["No population allele", "too few telomeric reads"]);
    });
    put_telomere(&d, 4000);
    with_reports(&d, |txt, html| {
        has(txt, &["Intratelomeric reads: 4,000"]);
        for report in [txt, html] {
            assert!(!report.contains("too few telomeric reads"));
        }
    });
}

#[test]
fn roh_reports_distinguish_an_upd_hint_from_a_relationship_band() {
    let d = Sample::new("report-roh-upd", "S");
    put_roh_with_sites(&d, UPD_RUNS);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &["chr15 carries 81.5 MB", "15q11-q13 (SNRPN, UBE3A)"],
            );
            assert!(!report.contains("Prader-Willi"));
            assert!(!report.contains("Angelman"));
        }
    });
    put_roh_with_sites(&d, UPD_DISTAL);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    "15q11-q13 (SNRPN, UBE3A)",
                    "F_ROH is the share of the autosomes",
                ],
            );
            for absent in ["The band is the row", "Prader-Willi", "Angelman"] {
                assert!(!report.contains(absent));
            }
        }
    });
    put_roh_with_sites(&d, UPD_OVER_HALF);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &["chr2 carries 25.1 MB, more than 50% of the length"],
            );
            assert!(!report.contains("MB, 50%"));
        }
    });
    put_roh_spread(&d);
    let f_roh = format!("F_ROH: {:.4}", 12.0 * SPREAD_LENGTH as f64 / ROH_DENOM);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &[
                    &f_roh,
                    "first cousins once removed (F 1/32)",
                    "The band is the row",
                ],
            );
            assert!(!report.contains("may come from one parent"));
        }
    });
}

#[test]
fn annotsv_reports_list_rule_downgrades_separately() {
    let d = Sample::new("report-annotsv-rules", "S");
    put_annotsv_checked(&d, false);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &[
                "step 22's consensus did not run",
                "support: not known (no step 22 consensus)",
            ],
        );
        for report in [txt, html] {
            assert!(!report.contains("Moved to class 3"));
        }
    });
    put_annotsv_checked(&d, true);
    with_reports(&d, |txt, html| {
        let note = og_render::notes::ANNOTSV_NOTE;
        assert!(note.contains("Step 05's rules"));
        has(txt, &[note]);
        has(
            html,
            &[
                "<h3>Moved to class 3 by step 05's rules</h3>",
                "<td>DEL 20:20000000-20450000</td><td>4</td><td>3</td><td>G9</td><td>ig_tr_locus</td>",
                "inside the IGFIX locus",
            ],
        );
    });
    fs::remove_file(d.path(ANNOTSV_CHECKED)).unwrap();
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            assert!(!report.contains("Class 4 and 5 (AnnotSV)"));
        }
    });
}

fn object(value: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    let Value::Object(fields) = value else {
        panic!("expected a summary object")
    };
    fields
}

#[test]
fn construct_hits_past_the_first_fifty_still_get_a_report_line() {
    let d = Sample::new("report-construct-past-limit", "S");
    put_construct(&d);
    let mut summary = d.collect();
    let sections = object(object(&mut summary).get_mut("sections").unwrap());
    let clinvar = object(sections.get_mut("clinvar").unwrap());
    let data = object(clinvar.get_mut("data").unwrap());
    let Value::Array(hits) = data.get_mut("hits").unwrap() else {
        panic!("missing ClinVar hits")
    };
    let cftr = hits.iter().find(|hit| matches!(hit, Value::Object(fields) if fields["gene"] == Value::Str("CFTR".into()))).unwrap().clone();
    let app = hits.iter().find(|hit| matches!(hit, Value::Object(fields) if fields["gene"] == Value::Str("APP".into()))).unwrap().clone();
    *hits = (100..160)
        .map(|pos| {
            let mut hit = cftr.clone();
            object(&mut hit).insert("pos".into(), Value::Int(pos));
            hit
        })
        .chain([app])
        .collect();
    render_value(&d, &summary, |txt, html| {
        let past = "APP chr21:25897626 TC>GA het: known lab-construct site (APP Swedish K670N/M671L): clear (0 of 12 alt reads and 0 of 10 ref reads clipped at an exon edge)";
        has(
            txt,
            &[
                "Hits at a known lab-construct site past the first 50:",
                past,
            ],
        );
        has(
            html,
            &[&og_render::escape::e(&format!("past the first 50: {past}"))],
        );
    });
}

#[test]
fn repeat_count_with_no_confidence_interval_does_not_size_an_in_repeat_call() {
    let locus =
        loads(r#"{"repeat_count":"20/300","allele_sizing":"sized/in_repeat","ci":"."}"#).unwrap();
    assert_eq!(
        og_render::sections::str::count(&locus).unwrap(),
        "20/not sized"
    );
}

#[test]
fn prs_effect_is_withheld_when_a_score_does_not_apply_to_the_declared_sex() {
    use og_render::sections::prs::effect;
    let mut score = loads(r#"{"effect":"HR 1.53 (EUR, PPM000591)","evaluated_in":["EUR"],"percentile":"81.2","group":"EUR","sex_check":"not applicable to the declared sex"}"#).unwrap();
    assert!(effect(&score).unwrap().0.is_empty());
    let fields = object(&mut score);
    fields.insert("sex_check".into(), Value::Str(String::new()));
    assert_eq!(effect(&score).unwrap().0, PRS_EUR_EFFECT);
}

#[test]
fn both_reports_name_cnvpytor_filters_and_counts() {
    let d = Sample::new("report-cnvpytor", "S");
    put_cnv_calls(&d);
    render_value(&d, &d.collect_with_sex(Some("female")), |txt, html| {
        has(
            txt,
            &[
                "Significant (e-val < 0.01): 5",
                "Significant after filters: 1 (left out: 1 q0 0.5 or more",
                "on chrY of a declared female",
            ],
        );
        has(html, &["CNVpytor significant after filters"]);
    });
    with_reports(&d, |txt, _| {
        has(txt, &["no declared sex, so chrY calls are kept"])
    });
}

#[test]
fn both_reports_name_only_plausible_karyotypes() {
    let d = Sample::new("report-karyotypes", "S");
    for (x, y, name) in KARYOTYPES {
        put_karyotype(&d, x, y);
        with_reports(&d, |txt, html| {
            if let Some(name) = name {
                has(
                    txt,
                    &[
                        &format!("Sex chromosomes: {name} (indexcov"),
                        &format!("Sex chromosomes: {name} (somalier"),
                    ],
                );
            }
            for report in [txt, html] {
                if let Some(name) = name {
                    has(report, &["Sex chromosomes", name]);
                } else {
                    assert!(!report.contains("Sex chromosomes"), "{x}, {y}");
                }
            }
        });
    }
}

#[test]
fn construct_site_notes_show_read_counts_when_measured() {
    use og_render::sections::{cdna::cdna_text, clinvar::D};
    let mut hit = loads(r#"{"cdna":{"check":"not_checked","note":"fewer than 2 reads carry the alt allele","construct":"LRRK2 G2019S","alt_reads":1,"alt_clipped":0,"ref_reads":6,"ref_clipped":0}}"#).unwrap();
    assert!(
        cdna_text(&D::new(&hit, "hit"))
            .unwrap()
            .ends_with("(0 of 1 alt reads and 0 of 6 ref reads clipped at an exon edge)")
    );
    let mut no_bam = hit.clone();
    let cdna = object(object(&mut no_bam).get_mut("cdna").unwrap());
    cdna.insert("note".into(), Value::Str("no BAM".into()));
    cdna.insert("alt_reads".into(), Value::Null);
    assert_eq!(
        cdna_text(&D::new(&no_bam, "hit")).unwrap(),
        "known lab-construct site (LRRK2 G2019S): not checked: no BAM"
    );
    let cdna = object(object(&mut hit).get_mut("cdna").unwrap());
    cdna.insert("construct".into(), Value::Str(String::new()));
    assert_eq!(
        cdna_text(&D::new(&hit, "hit")).unwrap(),
        "not checked: fewer than 2 reads carry the alt allele"
    );
}

fn put_prs_without_evaluations(d: &Sample) {
    d.put(
        PRS,
        "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tMatched_Pct\tPercentile\tAncestry_Group\tChrX\tAmbiguous\tApplies_To\tSex_Check\tInput\nCoronary artery disease\tPGS000018\t0.5\t1\t2\t50\t81.2\tEUR\t0\t0\tall\tapplies\tgvcf\n",
    );
}

const ANNOTSV_CHECKED: &str = "annotsv/S_sv_class4_5_checked.tsv";

fn put_annotsv_checked(d: &Sample, moved: bool) {
    // The checked table uses the columns written by og sv postfilter.
    let header = "ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tACMG_class\tCaller\tGene_name\tExon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\tAnnotSV_class\tRule\tRule_detail\n";
    let last = if moved {
        "SV3\t20\t20000000\t20450000\tDEL\t3\tmanta\tG9\tno\t\t\t\tno\t4\tig_tr_locus\tinside the IGFIX locus\n"
    } else {
        "SV3\t20\t20000000\t20450000\tDEL\t4\tmanta\tG9\tno\t\t\t\tno\t4\t\t\n"
    };
    d.put(ANNOTSV_CHECKED, format!("{header}SV1\t20\t1000200\t1009800\tDEL\t4\tmanta\tFIXG1\tyes\texon edges\t0.05\tgnomAD-SV_v4.1_DEL_fixture\tyes\t4\t\t\n\
        SV2\t20\t10100999\t10102999\tDEL\t5\tmanta\tG1;G2;G3;G4;G5\tno\t\t\t\tno\t5\t\t\n{last}"));
}

#[test]
fn acmg_secondary_findings_are_visible_in_both_reports() {
    let d = Sample::new("report-acmg-secondary-findings", "S");
    put_acmg_secondary_findings(&d);
    with_reports(&d, |txt, html| {
        has(
            txt,
            &["1 ClinVar P/LP, 1 rare HIGH impact", "    BRCA2", "    TTN"],
        );
        let html = card(html, "Secondary-Findings Genes (ACMG SF v3.3)");
        has(html, &["<td>TTN</td>"]);
        assert!(!html.contains("<td>GENEX</td>"));
    });
}

#[test]
fn blood_group_rows_use_the_same_schema_as_the_reports() {
    let d = Sample::new("report-blood-groups", "S");
    put_blood_groups(&d);
    with_reports(&d, |txt, html| {
        for report in [txt, html] {
            has(
                report,
                &["ABO", "A1/O.01", "RhD", "RHD deleted/RHD deleted", "Kell"],
            );
        }
    });
}

#[test]
fn corrupt_vcf_leaves_clinvar_readable_in_both_reports() {
    let d = Sample::new("report-corrupt-vcf", "S");
    put_bash(&d);
    let name = "vcf/S.vcf.gz";
    for (_, bytes) in damaged_vcfs(&d, name) {
        d.put(name, bytes);
        with_reports(&d, |txt, html| {
            has(txt, &["Pathogenic/Likely Pathogenic hits: 3"]);
            has(html, &["badge-yellow\">3</span>"]);
            has(card(html, "Variant Calling"), &["badge-gray\">Unreadable<"]);
        });
    }
}
