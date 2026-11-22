//! The collector cases from tests/test_collect_summary.py. The section map
//! lives in the PR alongside the cases already covered by owning crates.

#[path = "port_tests_summary_support/port_tests_summary_fixtures.rs"]
mod fixtures;

use fixtures::*;
use pycompat::json::{Dump, Value, dumps, loads};
use std::fs;

fn at<'a>(value: &'a Value, path: &str) -> &'a Value {
    path.split('.').fold(value, |value, key| match value {
        Value::Object(entries) => entries
            .get(key)
            .unwrap_or_else(|| panic!("missing {key} in {path}")),
        Value::Array(items) => &items[key.parse::<usize>().unwrap()],
        _ => panic!("cannot read {key} in {path}: {value:?}"),
    })
}

fn section<'a>(summary: &'a Value, key: &str) -> &'a Value {
    at(summary, &format!("sections.{key}"))
}

fn data<'a>(summary: &'a Value, key: &str) -> &'a Value {
    at(section(summary, key), "data")
}

fn expect(value: &Value, path: &str, json: &str) {
    assert_eq!(at(value, path), &loads(json).unwrap(), "{path}");
}

#[test]
fn bash_and_nextflow_folders_read_the_same_kinds_of_results() {
    let bash = Sample::new("bash", "S");
    put_bash(&bash);
    let summary = bash.collect();
    let clinvar = data(&summary, "clinvar");
    expect(clinvar, "count", "3");
    expect(clinvar, "hits.0.gene", r#""GENEB""#);
    expect(clinvar, "hits.1.gene", r#""GENEA""#);
    expect(clinvar, "hits.2.gene", r#""GENEC""#);
    assert_eq!(
        at(clinvar, "by_stars"),
        &loads(r#"{"4":0,"3":1,"2":0,"1":1,"0":1}"#).unwrap()
    );
    expect(data(&summary, "coverage"), "mean_depth", "29.5");
    expect(data(&summary, "sex_check"), "inferred_sex", r#""female""#);
    assert_eq!(
        data(&summary, "variants"),
        &loads(r#"{"total":3,"pass":2,"snps":2,"indels":1}"#).unwrap()
    );
    expect(data(&summary, "prs"), "scores.0.pgs_id", r#""PGS000018""#);
    expect(data(&summary, "prs"), "scores.1.pgs_id", r#""PGS000014""#);
    expect(data(&summary, "cpic"), "non_normal", "1");
    expect(data(&summary, "cpic"), "unclassified", "1");
    expect(&summary, "cyp2d6.agree", "true");
    expect(data(&summary, "mito"), "pass", "5");
    expect(data(&summary, "mito"), "heteroplasmic", "3");
    expect(
        data(&summary, "cpsr"),
        "classification",
        r#"{"VUS":2,"Likely_Pathogenic":1}"#,
    );

    let nextflow = Sample::new("nextflow", "S");
    put_nextflow(&nextflow);
    let summary = nextflow.collect();
    expect(section(&summary, "roh"), "state", r#""ok""#);
    expect(
        data(&summary, "roh"),
        "autosomal_over_5mb",
        r#"[{"region":"chr2:1-6000001","mb":6.0}]"#,
    );
    expect(data(&summary, "coverage"), "mean_depth", "31.25");
    expect(
        data(&summary, "hla"),
        "loci.0.alleles",
        r#"["B*07:02:01","B*08:01:01"]"#,
    );
    expect(data(&summary, "pharmcat"), "version", r#""3.2.0""#);
}

#[test]
fn the_run_status_distinguishes_stale_current_failed_and_missing() {
    let stale = Sample::new("stale", "S");
    put_stale(&stale);
    let summary = stale.collect();
    expect(section(&summary, "slivar"), "state", r#""stale""#);
    assert!(
        matches!(at(section(&summary,"slivar"), "note"), Value::Str(s) if s.contains("earlier run"))
    );
    expect(section(&summary, "roh"), "state", r#""ok""#);
    expect(section(&summary, "clinvar"), "state", r#""ok""#);
    expect(&summary, "run.declared_sex", r#""male""#);
    assert!(string(section(&summary, "slivar"), "note").contains("skipped"));
    fs::remove_file(stale.path(RUN)).unwrap();
    let summary = stale.collect();
    let Value::Object(sections) = at(&summary, "sections") else {
        panic!("missing sections")
    };
    assert!(
        sections
            .values()
            .all(|section| string(section, "state") != "stale")
    );

    let failed = Sample::new("failed", "S");
    put_failed(&failed);
    let summary = failed.collect();
    expect(section(&summary, "telomere"), "state", r#""failed""#);
    assert!(string(section(&summary, "telomere"), "note").contains("step 10 failed"));
    expect(section(&summary, "roh"), "state", r#""missing""#);
    expect(&summary, "cyp2d6.calls.pypgx", r#""failed (step 32)""#);
}

#[test]
fn malformed_gzip_is_unreadable_and_plain_text_with_a_gz_suffix_is_readable() {
    let d = Sample::new("gzip", "S");
    put_bash(&d);
    let name = "vcf/S.vcf.gz";
    let good = fs::read(d.path(name)).unwrap();
    let mut corrupt = good.clone();
    let crc = corrupt.len() - 8;
    corrupt[crc] ^= 1;
    for bytes in [b"\x1f\x8b\x08\x00not gzip\n".as_slice(), &corrupt] {
        d.put(name, bytes);
        let summary = d.collect();
        expect(section(&summary, "variants"), "state", r#""unreadable""#);
        expect(data(&summary, "clinvar"), "count", "3");
    }
    d.put(
        name,
        format!("{VCF}chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"),
    );
    expect(data(&d.collect(), "variants"), "total", "1");
}

#[test]
fn clinvar_cdna_checks_match_the_exact_allele() {
    let d = Sample::new("cdna-join", "S");
    put_bash(&d);
    d.put("clinvar/S_clinvar_cdna_check.tsv", "chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\nchr1\t10\tC\tT\t9\t8\t11\t0\tflagged\tpossible cDNA/transgene contamination: check the reads\nchr1\t20\tT\tA\t9\t0\t0\t0\tclear\t\n");
    let summary = d.collect();
    let cv = data(&summary, "clinvar");
    expect(cv, "cdna_flagged", "1");
    expect(cv, "hits.1.cdna.check", r#""flagged""#);
    assert!(string(cv, "hits.1.cdna.note").starts_with("possible cDNA"));
    for i in [0, 2] {
        expect(cv, &format!("hits.{i}.cdna.check"), r#""not_checked""#);
        expect(
            cv,
            &format!("hits.{i}.cdna.note"),
            r#""not in the read-check table""#,
        );
    }
}

#[test]
fn indexed_vep_maps_frequencies_and_rsids_to_the_hit_allele() {
    let d = Sample::new("clinvar-vep", "S");
    let hits = put_clinvar_quality(&d);

    let summary = d.collect();
    let cv = data(&summary, "clinvar");
    expect(cv, "count", "3");
    expect(cv, "low_confidence", "1");
    expect(cv, "by_stars", r#"{"4":0,"3":1,"2":0,"1":2,"0":0}"#);
    expect(cv, "frequency.source", r#""vep/S_vep.vcf.gz""#);
    expect(cv, "frequency.note", "null");
    let Value::Array(found) = at(cv, "hits") else {
        panic!("missing ClinVar hits")
    };
    assert_eq!(found.len(), 4);
    for (gene, af, rsids, low, common, gq, dp) in [
        ("GENEA", "0.03", r#"["rs123"]"#, "true", "true", "12", "6"),
        (
            "GENEB",
            "0.0002",
            r#"["rs9","rs77"]"#,
            "false",
            "false",
            "50",
            "40",
        ),
        ("GENEC", "0.0001", "[]", "false", "false", "45", "35"),
        ("GENED", "null", "[]", "false", "false", "null", "null"),
    ] {
        let hit = row(cv, "hits", "gene", gene);
        expect(hit, "af", af);
        expect(hit, "rsids", rsids);
        expect(hit, "low_confidence", low);
        expect(hit, "common", common);
        expect(hit, "gq", gq);
        expect(hit, "dp", dp);
    }
    let common = row(cv, "hits", "gene", "GENEA");
    expect(common, "gnomad_exome_af", "0.03");
    expect(common, "gnomad_genome_af", "0.02");

    let mut hit_rows = Vec::new();
    og_clinvar::hits::hits(hits.as_bytes(), &mut hit_rows).unwrap();
    let hit_rows = String::from_utf8(hit_rows).unwrap();
    let mut columns: Vec<_> = hit_rows
        .lines()
        .map(|line| {
            let cells: Vec<_> = line.split('\t').collect();
            (cells[6], cells[10], cells[11], cells[12])
        })
        .collect();
    columns.sort_unstable();
    assert_eq!(
        columns,
        [
            ("GENEA", "6", "12", "yes"),
            ("GENEB", "40", "50", "no"),
            ("GENEC", "35", "45", "no"),
            ("GENED", ".", ".", "no"),
        ]
    );

    let history = "history/20200101T000000/summary.json";
    let mut old = summary.clone();
    let Value::Object(root) = &mut old else {
        unreachable!()
    };
    let Some(Value::Object(sections)) = root.get_mut("sections") else {
        unreachable!()
    };
    let Some(Value::Object(clinvar)) = sections.get_mut("clinvar") else {
        unreachable!()
    };
    let Some(Value::Object(clinvar_data)) = clinvar.get_mut("data") else {
        unreachable!()
    };
    let Some(Value::Array(rows)) = clinvar_data.get_mut("hits") else {
        unreachable!()
    };
    for hit in rows.iter_mut() {
        let Value::Object(row) = hit else {
            unreachable!()
        };
        assert!(row.shift_remove("rsids").is_some());
    }
    d.put(history, dumps(&old, &Dump::default()).unwrap());
    expect(&d.collect(), "new_since.sections", "[]");

    let reclassified = hits.replacen("CLNSIG=Pathogenic", "CLNSIG=Likely_pathogenic", 1);
    d.put(CLINVAR, reclassified);
    expect(&d.collect(), "new_since.sections", r#"["clinvar"]"#);
    fs::remove_file(d.path(history)).unwrap();

    let zero = hits.replace("GQ=30,DP=10", "GQ=0,DP=0");
    d.put(CLINVAR, &zero);
    let mut zero_rows = Vec::new();
    og_clinvar::hits::hits(zero.as_bytes(), &mut zero_rows).unwrap();
    assert!(!String::from_utf8(zero_rows).unwrap().contains("\tyes"));
    expect(data(&d.collect(), "clinvar"), "count", "4");
    expect(data(&d.collect(), "clinvar"), "low_confidence", "0");
    fs::remove_file(d.path(&format!("{VEP}.tbi"))).unwrap();
    let no_index = d.collect();
    let cv = data(&no_index, "clinvar");
    assert!(matches!(at(cv, "frequency.note"), Value::Str(s) if s.contains("no .tbi")));
    fs::remove_dir_all(d.path("vep")).unwrap();
    let no_vep = d.collect();
    for cv in [cv, data(&no_vep, "clinvar")] {
        for i in 0..4 {
            expect(cv, &format!("hits.{i}.af"), "null");
            expect(cv, &format!("hits.{i}.rsids"), "[]");
        }
    }
}

#[test]
fn cpsr_compares_both_class_columns_and_joins_its_read_checks() {
    let d = Sample::new("cpsr", "S");
    put_cpsr(&d);
    let summary = d.collect();
    let cpsr = data(&summary, "cpsr");
    expect(
        cpsr,
        "classification",
        r#"{"Pathogenic":1,"Likely Pathogenic":2,"Benign":1,"VUS":1}"#,
    );
    expect(cpsr, "final_plp_own_benign", "2");
    expect(cpsr, "low_confidence", "1");
    expect(cpsr, "cdna_flagged", "1");
    for (i, gene, check) in [
        (0, "BRCA2", "clear"),
        (1, "BRCA1", "not_checked"),
        (2, "MSH2", "flagged"),
        (3, "APC", "not_checked"),
    ] {
        expect(cpsr, &format!("plp.{i}.gene"), &format!("\"{gene}\""));
        expect(
            cpsr,
            &format!("plp.{i}.cdna.check"),
            &format!("\"{check}\""),
        );
    }
    expect(cpsr, "plp.1.low_confidence", "true");
    expect(cpsr, "disagreements.0.acmg_codes", r#""BA1""#);
    expect(cpsr, "disagreements.1.acmg_codes", r#""BS1""#);
    expect(cpsr, "disagreements.2.acmg_codes", r#""PVS1,PM2""#);
    let disagreements = array(cpsr, "disagreements");
    assert_eq!(disagreements.len(), 3);
    assert_eq!(
        disagreements
            .iter()
            .map(|item| (
                string(item, "gene"),
                string(item, "clinvar"),
                string(item, "cpsr")
            ))
            .collect::<Vec<_>>(),
        [
            ("BRCA2", "Pathogenic", "Benign"),
            ("BRCA1", "Likely Pathogenic", "Likely Benign"),
            ("APC", "Benign", "Pathogenic"),
        ]
    );
    expect(&summary, "databases.cpsr_clinvar_release", r#""20250309""#);

    d.put(
        CLINVAR,
        VCF.replacen("#CHROM", "##PGP_CALL_FLOOR=GQ=0,DP=0\n#CHROM", 1),
    );
    expect(data(&d.collect(), "cpsr"), "low_confidence", "0");

    d.gzip("cpsr/S.cpsr.grch38.classification.tsv.gz", "VAR_ID\tSYMBOL\tCLASSIFICATION\n13_32340300_G_A\tBRCA2\tPathogenic\n5_112800000_T_C\tAPC\tBenign\n");
    let summary = d.collect();
    expect(data(&summary, "cpsr"), "disagreements", "[]");
    expect(data(&summary, "cpsr"), "own_class_column", "null");
}

#[test]
fn read_backed_phase_is_kept_in_both_folder_layouts() {
    for (name, nextflow, counts) in [
        ("pairs-bash", false, r#"{"trans":1,"cis":1,"unknown":1}"#),
        ("pairs-nextflow", true, r#"{"trans":1,"cis":0,"unknown":0}"#),
    ] {
        let d = Sample::new(name, "S");
        put_comphet(&d, nextflow);
        let summary = d.collect();
        let slivar = data(&summary, "slivar");
        expect(slivar, "compound_het_phase", counts);
        expect(
            slivar,
            "compound_het_trans",
            r#"[{"gene":"G2","variants":"chr1:50 A>G + chr1:90 C>T","reads":"19 of 20"}]"#,
        );
        if nextflow {
            expect(slivar, "prioritized", "1");
        }
    }
}

#[test]
fn couple_rows_follow_the_sample_and_reject_another_couple() {
    let d = Sample::new("couple", "B");
    let empty = d.collect();
    expect(section(&empty, "couple"), "state", r#""missing""#);
    assert!(
        matches!(at(&empty, "not_run"), Value::Array(items) if items.iter().all(|v| !matches!(v, Value::Str(s) if s.contains("Couple"))))
    );
    put_couple(&d);
    let summary = d.collect();
    let couple = data(&summary, "couple");
    expect(section(&summary, "couple"), "step", r#""39""#);
    expect(section(&summary, "couple"), "state", r#""ok""#);
    assert_eq!(array(couple, "shared").len(), 1);
    expect(couple, "partner", r#""A""#);
    expect(couple, "shared.0.gene", r#""NAGA""#);
    expect(couple, "shared.0.own", r#""b-snv""#);
    expect(couple, "shared.0.partner_findings", r#""a-snv""#);
    expect(couple, "x_linked.0.gene", r#""FMR1""#);
    expect(couple, "one_partner", "2");
    let body = fs::read_to_string(d.path(COUPLE)).unwrap();
    for invalid in [
        body.replace("\tB\t", "\tC\t"),
        "gene\tresult\nNAGA\tshared\n".into(),
    ] {
        d.put(COUPLE, invalid);
        expect(section(&d.collect(), "couple"), "state", r#""unreadable""#);
    }
    d.put(COUPLE, COUPLE_HEADER);
    let empty_table = d.collect();
    expect(section(&empty_table, "couple"), "state", r#""ok""#);
    expect(data(&empty_table, "couple"), "partner", "null");
    expect(data(&empty_table, "couple"), "shared", "[]");
}

#[test]
fn export_counts_are_required_and_join_the_run_status() {
    let d = Sample::new("export", "S");
    put_export(&d);
    let summary = d.collect();
    expect(section(&summary, "consumer_export"), "step", r#""40""#);
    expect(section(&summary, "consumer_export"), "state", r#""ok""#);
    let export = data(&summary, "consumer_export");
    expect(export, "sites", "10");
    expect(export, "called", "7");
    expect(export, "no_call_reasons", r#"{"no-call in the gVCF":3}"#);
    for invalid in [
        "key\tvalue\nsites\t10\ncalled\t\nno_call\t3\n",
        "key\tvalue\nsites\t10\ncalled\t7\nno_call\t3\nno_call: filtered call\tthree\n",
    ] {
        d.put(EXPORT, invalid);
        expect(
            section(&d.collect(), "consumer_export"),
            "state",
            r#""unreadable""#,
        );
    }
    put_export(&d);
    d.old(EXPORT);
    d.put(
        RUN,
        format!("meta\tstarted_epoch\t{STARTED}\nstep\t40\tskipped (opt-in)\n"),
    );
    expect(
        section(&d.collect(), "consumer_export"),
        "state",
        r#""stale""#,
    );
}

#[test]
fn mosaic_estimates_share_the_qc_table_and_refuse_an_unknown_call() {
    let d = Sample::new("mosaic", "S");
    d.put(
        QC,
        "key\tvalue\nsample\tS\ninferred_sex\tmale\nsex_check\tok\n",
    );
    expect(
        section(&d.collect(), "mosaic_sex_loss"),
        "state",
        r#""missing""#,
    );
    put_mosaic(&d);
    let summary = d.collect();
    let mosaic = section(&summary, "mosaic_sex_loss");
    expect(mosaic, "step", r#""33""#);
    expect(mosaic, "source", r#""qc/S_sample_qc.tsv""#);
    expect(mosaic, "data.loss_depth", "0.298");
    expect(mosaic, "data.loss_af", "null");
    expect(mosaic, "data.x_het_sites", "null");
    expect(mosaic, "data.call", r#""detected""#);
    expect(mosaic, "data.chrom", r#""chrY""#);
    expect(data(&summary, "sample_qc"), "sex_check", r#""ok""#);
    let rows = fs::read_to_string(d.path(QC)).unwrap();
    d.put(QC, rows.replace("0.298", "inf").replace("0.303", "-inf"));
    let nonfinite = d.collect();
    expect(data(&nonfinite, "mosaic_sex_loss"), "loss_depth", "null");
    expect(
        data(&nonfinite, "mosaic_sex_loss"),
        "loss_depth_high",
        "null",
    );
    d.put(
        QC,
        rows.replace("mosaic_call\tdetected", "mosaic_call\tmaybe"),
    );
    expect(
        section(&d.collect(), "mosaic_sex_loss"),
        "state",
        r#""unreadable""#,
    );
}

#[test]
fn unusable_pgx_calls_and_stale_consensus_do_not_claim_agreement() {
    let d = Sample::new("unusable-cyp", "S");
    put_cyp_unusable(&d);
    let summary = d.collect();
    expect(data(&summary, "pypgx"), "genes_called", "1");
    expect(&summary, "cyp2d6.agree", "null");
    expect(
        &summary,
        "cyp2d6.calls.Cyrius",
        r#""*5/*5 (not usable: Filter CYP2D6_depth_unreliable)""#,
    );
    expect(&summary, "cyp2d6.consensus.result", r#""indeterminate""#);
    expect(&summary, "cyp2d6.consensus.passed_to_pharmcat", "false");
    let verdict = fs::read(d.path(CONSENSUS)).unwrap();
    fs::remove_file(d.path(CONSENSUS)).unwrap();
    expect(&d.collect(), "cyp2d6.consensus", "null");
    d.put(CONSENSUS, verdict);
    d.old(CONSENSUS);
    expect(&d.collect(), "cyp2d6.consensus", "null");
}

#[test]
fn prs_keeps_the_group_where_the_effect_was_evaluated() {
    for (group, effect) in [("EUR", PRS_EUR_EFFECT), ("EAS", "")] {
        let d = Sample::new(group, "S");
        put_prs_evaluation(&d, group, effect);
        let summary = d.collect();
        let row = at(data(&summary, "prs"), "scores.0");
        expect(row, "evaluated_in", r#"["EUR","AFR","AMR"]"#);
        expect(row, "effect", &format!("\"{effect}\""));
    }
}

#[test]
fn lineage_reads_both_y_placements_and_the_contamination_check() {
    let d = Sample::new("lineage", "S");
    put_lineage(&d);
    let summary = d.collect();
    expect(
        data(&summary, "haplogroup"),
        "contamination_status",
        r#""YES""#,
    );
    expect(
        data(&summary, "haplogroup"),
        "contamination_level",
        r#""0.12""#,
    );
    expect(data(&summary, "y_haplogroup"), "haplogroup", r#""R-M269""#);
    let y = data(&summary, "y_tree");
    expect(y, "haplogroup", r#""R-FT12345""#);
    expect(y, "tree", r#""YFull YTree 14.06.0""#);
    expect(y, "level", r#""24""#);
    expect(y, "path.0.Branch", r#""A0-T""#);
    expect(y, "path.1.Branch", r#""R-FT12345""#);
    put_unplaced_y(&d);
    let summary = d.collect();
    expect(
        data(&summary, "y_haplogroup"),
        "haplogroup",
        r#""insufficient markers""#,
    );
    expect(section(&summary, "y_tree"), "state", r#""ok""#);
    let y = data(&summary, "y_tree");
    expect(y, "haplogroup", r#""no branch has a derived call""#);
    expect(y, "path", "[]");
    d.put("y_haplogroup/S_y_tree_path.tsv", "Level\tSNPs\n");
    expect(section(&d.collect(), "y_tree"), "state", r#""unreadable""#);
}

#[test]
fn stranger_keeps_ungraded_loci_and_fresh_motifs_separate() {
    let d = Sample::new("repeat-followups", "S");
    put_ungraded_repeat(&d);
    let summary = d.collect();
    let expansions = data(&summary, "expansions");
    expect(
        expansions,
        "not_graded.MUC1",
        r#""the reference repeat is 1200 bases, longer than 300: short reads cannot size it""#,
    );
    expect(expansions, "no_status", "1");
    let muc1 = row(expansions, "loci", "locus", "MUC1");
    assert_eq!(at(muc1, "not_graded"), at(expansions, "not_graded.MUC1"));
    expect(
        expansions,
        "flagged.0.motif.verdict",
        r#""pathogenic_motif""#,
    );
    d.old("expansion_hunter/S_str_motifs.tsv");
    let old = d.collect();
    let hit = at(data(&old, "expansions"), "flagged.0");
    assert!(matches!(hit, Value::Object(row) if !row.contains_key("motif")));
    fs::remove_file(d.path(STR_CATALOG)).unwrap();
    let no_catalog = d.collect();
    expect(data(&no_catalog, "expansions"), "not_graded", "{}");
    let expansions = data(&no_catalog, "expansions");
    expect(expansions, "no_status", "2");
    assert!(
        array(expansions, "loci").iter().all(
            |item| matches!(item, Value::Object(fields) if !fields.contains_key("not_graded"))
        )
    );
}

#[test]
fn the_partial_call_construct_scan_and_always_report_rows_reach_the_summary() {
    let d = Sample::new("construct", "S");
    put_construct(&d);
    let summary = d.collect();
    let cv = data(&summary, "clinvar");
    expect(cv, "hits.1.gene", r#""CFTR""#);
    expect(
        cv,
        "hits.1.genotype",
        r#""one copy, other copy not called ref (1/.)""#,
    );
    expect(cv, "hits.0.cdna.construct", r#""APP Swedish K670N/M671L""#);
    expect(cv, "hits.0.cdna.alt_reads", "12");
    expect(cv, "hits.0.cdna.ref_reads", "10");
    expect(cv, "hits.0.cdna.check", r#""clear""#);
    expect(cv, "construct_sites", "1");
    expect(
        data(&summary, "carrier"),
        "carriers.0.genotype",
        r#""one copy, other copy not called ref (1/.)""#,
    );

    let mut hits = Vec::new();
    og_clinvar::hits::hits(
        std::io::BufReader::new(fs::File::open(d.path(CLINVAR)).unwrap()),
        &mut hits,
    )
    .unwrap();
    let hits = String::from_utf8(hits).unwrap();
    assert_eq!(
        hits.lines()
            .map(|row| row.split('\t').nth(5).unwrap())
            .collect::<Vec<_>>(),
        ["one copy, other copy not called ref (1/.)", "het"]
    );
}

#[test]
fn paralog_counts_keep_the_exon_side_quality_and_not_assessed_state() {
    let d = Sample::new("parascopy", "S");
    put_parascopy(&d);
    let summary = d.collect();
    let p = data(&summary, "parascopy");
    expect(section(&summary, "parascopy"), "state", r#""ok""#);
    assert!(string(p, "reading").starts_with("two SMN1 copies: not a carrier"));
    assert!(string(p, "reading").contains("2+0"));
    expect(section(&summary, "parascopy"), "step", r#""35""#);
    expect(p, "smn1.cn", "2");
    expect(p, "smn2.cn", "2");
    expect(p, "agcn", r#""4""#);
    expect(p, "agcn_qual", r#""60""#);
    expect(p, "hba.alpha_globin_genes", r#""4""#);
    expect(p, "smn1_marker.status", r#""present""#);

    d.put(SMN, format!("{SMN_HEADER}{SMN_CARRIER}"));
    let summary = d.collect();
    expect(data(&summary, "parascopy"), "smn1.cn", "1");
    expect(data(&summary, "parascopy"), "smn2.cn", "2");
    expect(
        data(&summary, "parascopy"),
        "reading",
        r#""one SMN1 copy: the pattern of an SMA carrier""#,
    );
    d.put(SMN, format!("{SMN_HEADER}chr5\t70950920\t70953320\tSMN1\tPASS\t3\t35\tPASS\t1,2\t15,16\tchr5:70075496-70077900:+\n"));
    let weak = d.collect();
    let p = data(&weak, "parascopy");
    expect(p, "smn1.cn", "null");
    expect(
        p,
        "smn1.text",
        r#""1 (quality 15, filter PASS), not resolved""#,
    );
    assert!(string(p, "reading").starts_with("SMN1 not resolved"));

    d.put(SMN, format!("{SMN_HEADER}{SMN_ELSEWHERE}"));
    expect(
        data(&d.collect(), "parascopy"),
        "reading",
        r#""no row of step 35's table covers SMN1 exon 7: no carrier reading""#,
    );
    d.put(SMN, format!("{SMN_HEADER}chr5\n{SMN_EXON_7}"));
    expect(data(&d.collect(), "parascopy"), "smn1.cn", "2");

    d.put(SMN, "chrom\tstart\tend\tagCN\nchr5\t1\t2\t4\n");
    expect(
        section(&d.collect(), "parascopy"),
        "state",
        r#""unreadable""#,
    );
}

#[test]
fn clinical_mane_and_call_quality_keep_weak_and_unknown_rows_out_of_the_rare_count() {
    let d = Sample::new("clinical-mane", "S");
    put_clinical_mane(&d, false);
    let summary = d.collect();
    let cl = data(&summary, "clinical");
    let sf = at(cl, "acmg_sf");
    expect(cl, "frequency_unknown", "1");
    expect(cl, "weak_calls", "1");
    expect(data(&summary, "slivar"), "frequency_unknown", "1");
    expect(sf, "high_impact_rare", "2");
    expect(sf, "high_impact_weak", "2");
    expect(sf, "high_impact_frequency_unknown", "1");

    let Value::Array(rows) = at(sf, "high_impact") else {
        panic!("missing high impact rows")
    };
    let find = |gene: &str| {
        rows.iter()
            .find(|row| matches!(at(row, "gene"), Value::Str(g) if g == gene))
    };
    assert!(find("TTN").is_none() && find("KCNQ1").is_none());
    expect(find("SCN5A").unwrap(), "transcript", r#""CANONICAL""#);
    expect(find("SCN5A").unwrap(), "consequence", r#""stop_gained""#);
    expect(find("BRCA2").unwrap(), "weak_call", "true");
    expect(find("LDLR").unwrap(), "weak_call", "true");
    expect(find("MYBPC3").unwrap(), "weak_call", "false");
    expect(find("MYH7").unwrap(), "frequency", r#""unknown""#);
    expect(find("MYH7").unwrap(), "frequency_unknown", "true");
    assert!(string(find("BRCA2").unwrap(), "call_note").contains("GQ 12 < 30"));
    let ldlr = string(find("LDLR").unwrap(), "call_note");
    assert!(ldlr.contains("GQ 11") && ldlr.contains("0.29 of the depth on a homozygous call"));

    put_clinical_mane(&d, true);
    let summary = d.collect();
    let sf = at(data(&summary, "clinical"), "acmg_sf");
    expect(sf, "high_impact_rare", "4");
    expect(sf, "high_impact_weak", "0");
}

#[test]
fn annotsv_postfilter_preserves_the_phase_used_by_two_hits() {
    const INPUT: &str = "annotsv/S_sv_class4_5.tsv";
    const OUTPUT: &str = "annotsv/S_sv_class4_5_checked.tsv";
    let d = Sample::new("annotsv-phase", "S");
    d.put(CLINVAR, format!("{VCF}chr1\t45332000\t.\tC\tT\t50\tPASS\tGENEINFO=MUTYH:4595;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT:GQ:DP:PS\t1|0:50:30:45330000\n"));
    for (gt, verdict, reason) in [
        ("0|1", "biallelic", "in trans"),
        ("0/0", "", ""),
        (".|1", "possibly biallelic", "phase unknown"),
    ] {
        d.put(INPUT, format!(
            "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tID\tINFO\tFORMAT\tS\tAnnotation_mode\tGene_name\tB_loss_AFmax\tACMG_class\tAnnotSV_ranking_score\n\
             1_45335001_45337000_DEL_1\t1\t45335000\t45337000\t-2000\tDEL\tSV1\tEND=45337000;SVTYPE=DEL;SVCALLER=manta\tGT:PS\t{gt}:45330000\tfull\tMUTYH\t\t4\t0.95\n"
        ));
        assert_eq!(
            og_steps::sv::postfilter::cli(&[
                "--class45".into(),
                d.path(INPUT).display().to_string(),
                "--out".into(),
                d.path(OUTPUT).display().to_string(),
            ]),
            0,
        );
        let checked = fs::read_to_string(d.path(OUTPUT)).unwrap();
        assert!(checked.lines().next().unwrap().ends_with("\tFORMAT\tS"));
        assert!(checked.contains(&format!("GT:PS\t{gt}:45330000")));
        let summary = d.collect();
        let Value::Array(genes) = at(data(&summary, "carrier"), "two_hits.genes") else {
            panic!("missing two-hit gene list");
        };
        let pair = genes
            .iter()
            .find(|row| matches!(at(row, "gene"), Value::Str(g) if g == "MUTYH"));
        if verdict.is_empty() {
            assert!(pair.is_none());
        } else {
            let pair = pair.expect("the MUTYH deletion and ClinVar hit should pair");
            assert_eq!(at(pair, "verdict"), &Value::Str(verdict.into()));
            assert!(matches!(at(pair, "why"), Value::Str(why) if why.contains(reason)));
        }
    }
}

#[test]
fn blood_groups_require_the_phenotype_column_and_follow_step_41() {
    let d = Sample::new("blood-groups", "S");
    let summary = d.collect();
    assert!(
        matches!(at(&summary, "not_run"), Value::Array(v) if v.iter().any(|x| matches!(x, Value::Str(s) if s == "Blood groups (predicted)")))
    );
    put_blood_groups(&d);
    let summary = d.collect();
    expect(section(&summary, "bloodtype"), "state", r#""ok""#);
    expect(section(&summary, "bloodtype"), "step", r#""41""#);
    assert!(og_summary::collect::model::RUN_ALL_STEPS.contains(&"41"));
    expect(data(&summary, "bloodtype"), "abo", r#""A1""#);
    assert_eq!(array(data(&summary, "bloodtype"), "systems").len(), 3);
    expect(data(&summary, "bloodtype"), "rhd", r#""D-""#);
    expect(data(&summary, "bloodtype"), "systems.2.system", r#""Kell""#);
    d.put(BLOODTYPE, "system\tantigens\nABO\tA\n");
    let summary = d.collect();
    expect(section(&summary, "bloodtype"), "state", r#""unreadable""#);
    assert!(
        matches!(at(section(&summary, "bloodtype"), "note"), Value::Str(s) if s.contains("phenotype"))
    );
}

#[test]
fn roh_excludes_sex_chromosomes_and_masked_bands() {
    let d = Sample::new("masked-roh", "S");
    put_masked_roh(&d);
    let summary = d.collect();
    let roh = data(&summary, "roh");
    expect(roh, "largest_mb", "7.0");
    expect(roh, "total_mb", "10.5");
    expect(roh, "segments", "2");
    expect(
        roh,
        "autosomal_over_5mb",
        r#"[{"region":"chr2:10000001-17000000","mb":7.0}]"#,
    );
    expect(roh, "sex_chromosome_segments", "1");
    expect(roh, "masked_mb", "46.5");
    for path in ["scripts/11-roh-analysis.sh", "modules/local/roh/main.nf"] {
        let source = fs::read_to_string(std::path::Path::new(REPO).join(path)).unwrap();
        assert!(!source.contains("MASKED_BANDS="), "{path}");
        assert!(source.contains("og roh summary"), "{path}");
    }
    let args = [
        "--sample".into(),
        "S".into(),
        "--roh".into(),
        d.path(ROH).display().to_string(),
    ];
    let mut out = Vec::new();
    let mut err = Vec::new();
    assert_eq!(
        og_summary::collect::group_c::roh::cli(&args, &mut out, &mut err),
        0,
        "{}",
        String::from_utf8_lossy(&err)
    );
    assert_eq!(
        String::from_utf8(out)
            .unwrap()
            .lines()
            .skip(3)
            .collect::<Vec<_>>(),
        ["chr2\t10000001\t17000000\t7000000\t7.0"]
    );
    d.put(ROH, "RG\tS\tchr2\t10000001\t13000000\tinf\t900\t50\n");
    expect(section(&d.collect(), "roh"), "state", r#""unreadable""#);
}

#[test]
fn telomere_counts_require_enough_reads_at_the_observed_depth() {
    let d = Sample::new("telomere", "S");
    for (count, low) in [(40, "true"), (4000, "false")] {
        put_telomere(&d, count);
        let summary = d.collect();
        let telomere = data(&summary, "telomere");
        expect(telomere, "intratelomeric_reads", &count.to_string());
        expect(telomere, "low_reads", low);
        expect(telomere, "depth", "30.0");
        expect(telomere, "read_floor", "2000");
    }
    d.put(
        TELOMERE,
        "PID\tsample\tintratelomeric_reads\ttel_content\nS\ttumor\t61873\t334\n",
    );
    let summary = d.collect();
    let telomere = data(&summary, "telomere");
    expect(telomere, "intratelomeric_reads", "61873");
    assert!(matches!(at(telomere, "read_check"), Value::Str(s) if s.contains("not checked")));
    assert!(matches!(telomere, Value::Object(o) if !o.contains_key("low_reads")));
}

#[test]
fn the_readers_agree_on_clinvar_stars() {
    let statuses = [
        ("practice_guideline", 4),
        ("reviewed_by_expert_panel", 3),
        ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
        ("criteria_provided,_single_submitter", 1),
        ("criteria_provided,_conflicting_interpretations", 1),
        ("criteria_provided,_conflicting_classifications", 1),
        ("no_assertion_criteria_provided", 0),
        ("no_classification_provided", 0),
        ("something_new", 0),
        ("", 0),
    ];
    let d = Sample::new("clinvar-stars", "S");
    let mut input = VCF.to_owned();
    for (i, (status, _)) in statuses.iter().enumerate() {
        use std::fmt::Write;
        writeln!(input, "chr1\t{}\t.\tA\tG\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT={status}\tGT\t0/1", i + 1).unwrap();
    }
    d.put(CLINVAR, &input);
    let mut hits = Vec::new();
    og_clinvar::hits::hits(input.as_bytes(), &mut hits).unwrap();
    let got = String::from_utf8(hits).unwrap();
    assert_eq!(
        got.lines()
            .map(|line| line.split('\t').next().unwrap().parse::<i128>().unwrap())
            .collect::<Vec<_>>(),
        statuses.map(|(_, stars)| stars),
    );
    let summary = d.collect();
    let hits = array(data(&summary, "clinvar"), "hits");
    let mut actual: Vec<_> = hits
        .iter()
        .map(|hit| {
            let (Value::Int(pos), Value::Int(stars)) = (at(hit, "pos"), at(hit, "stars")) else {
                panic!("position and stars must be integers");
            };
            (*pos, *stars)
        })
        .collect();
    actual.sort_unstable();
    let expected: Vec<_> = statuses
        .iter()
        .enumerate()
        .map(|(i, (_, stars))| (i as i128 + 1, *stars))
        .collect();
    assert_eq!(actual, expected);
}

#[test]
fn clinvar_combined_allele_notes_are_exact_and_keep_star_order() {
    use std::collections::BTreeMap;
    let d = Sample::new("clinvar-combined-alleles", "S");
    for (gt, suffix) in [("0/1", COMBINED_UNPHASED), ("0|1", ""), ("1/1", "")] {
        let rows = combined_hits(gt);
        d.put(CLINVAR, &rows);
        let mut hit_rows = Vec::new();
        og_clinvar::hits::hits(rows.as_bytes(), &mut hit_rows).unwrap();
        let hit_rows = String::from_utf8(hit_rows).unwrap();
        let cli_notes: BTreeMap<_, _> = hit_rows
            .lines()
            .map(|line| {
                let cells: Vec<_> = line.split('\t').collect();
                (cells[2].parse::<i128>().unwrap(), cells[9].to_owned())
            })
            .collect();
        let summary = d.collect();
        let hits = array(data(&summary, "clinvar"), "hits");
        expect(&hits[0], "pos", "30000038");
        let notes: BTreeMap<_, _> = hits
            .iter()
            .map(|hit| {
                let Value::Int(pos) = at(hit, "pos") else {
                    panic!("invalid position")
                };
                (*pos, string(hit, "combined_allele").to_owned())
            })
            .collect();
        let part = format!("part of a combined allele: chr22:30000038 GCCAT>CATGGG{suffix}");
        let want = BTreeMap::from([
            (30000037, part.clone()),
            (
                30000038,
                format!("combined allele of chr22:30000037 AGC>A + chr22:30000042 T>TGGG{suffix}"),
            ),
            (30000042, part),
            (30000210, ".".into()),
            (
                30000300,
                "possibly one allele with chr22:30000305 C>T: check ClinVar for a combined record"
                    .into(),
            ),
            (
                30000305,
                "possibly one allele with chr22:30000300 A>G: check ClinVar for a combined record"
                    .into(),
            ),
        ]);
        assert_eq!(notes, want, "{gt}");
        assert_eq!(cli_notes, want, "{gt}");
    }
}

#[test]
fn hla_concordance_and_kir_keep_the_comparison_and_no_read_reasons() {
    let bash = Sample::new("hla-qc", "S");
    put_hla_qc(&bash, false);
    let expected = fs::read_to_string(bash.path(HLA_CONCORDANCE)).unwrap();
    fs::remove_file(bash.path(HLA_CONCORDANCE)).unwrap();
    og_steps::pgx_rules::hla_concordance(
        bash.path(HLA_GENOTYPES).to_str().unwrap(),
        bash.path(OPTITYPE).to_str().unwrap(),
        bash.path(HLA_CONCORDANCE).to_str().unwrap(),
    )
    .unwrap();
    let written = fs::read_to_string(bash.path(HLA_CONCORDANCE)).unwrap();
    assert_eq!(
        written.lines().collect::<Vec<_>>(),
        expected.lines().collect::<Vec<_>>()
    );
    let summary = bash.collect();
    let hla = data(&summary, "hla");
    for (index, gene, agree) in [
        (0, "HLA-A", "yes"),
        (1, "HLA-B", "no"),
        (2, "HLA-C", "not compared"),
    ] {
        expect(
            hla,
            &format!("optitype.genes.{index}.gene"),
            &format!("\"{gene}\""),
        );
        expect(
            hla,
            &format!("optitype.genes.{index}.agree"),
            &format!("\"{agree}\""),
        );
    }
    expect(hla, "kir.typed.0.gene", r#""KIR2DL1""#);
    expect(hla, "kir.typed.1.gene", r#""KIR3DL1""#);
    expect(hla, "kir.typed.0.quality", r#"["60","35"]"#);
    expect(hla, "kir.typed.1.alleles", r#"["KIR3DL1*0150201"]"#);
    expect(hla, "kir.not_called", r#"["KIR2DS1"]"#);

    let nextflow = Sample::new("hla-qc-nextflow", "S");
    put_hla_qc(&nextflow, true);
    let summary = nextflow.collect();
    let hla = data(&summary, "hla");
    expect(hla, "optitype.genes", "[]");
    expect(
        hla,
        "optitype.reason",
        r#""OptiType not run: no read pairs at HLA-A, HLA-B or HLA-C in this BAM""#,
    );
    expect(hla, "kir.typed", "[]");
    expect(
        hla,
        "kir.reason",
        r#""KIR not typed: T1K found too few reads at the KIR genes in this BAM""#,
    );
}

#[test]
fn pypgx_copy_number_and_probable_cyp2d6_keep_their_evidence() {
    let d = Sample::new("pgx-evidence", "S");
    put_pgx_evidence(&d, "ok");
    let summary = d.collect();
    let pgx = data(&summary, "pypgx");
    expect(
        pgx,
        "cnv_call",
        r#"{"CYP2D6":"Normal","CYP2A6":"Normal","GSTM1":"WholeDel2"}"#,
    );
    expect(pgx, "not_callable", r#"["GSTT1"]"#);
    expect(pgx, "genes_called", "4");
    assert!(
        matches!(data(&summary, "cnvpytor"), Value::Object(fields) if !fields.contains_key("ugt2b17"))
    );
    expect(
        pgx,
        "low_gq",
        r#"[{"gene":"CYP2C19","diplotype":"*1/*2","min_gq":12,"weakest":"chr10:94781859 G>A GQ 12 AD 26,3"}]"#,
    );
    expect(
        &summary,
        "cyp2d6.probable",
        r#"{"diplotype":"*1/*2","phenotype":"Normal Metabolizer","label":"probable (one caller: pypgx)"}"#,
    );
    put_pgx_evidence(&d, "unreliable");
    expect(&d.collect(), "cyp2d6.probable", "null");
    put_pgx_evidence(&d, "ok");
    put_cyrius(&d, "*1/*2", "PASS");
    expect(&d.collect(), "cyp2d6.probable", "null");
    d.put(
        CYRIUS,
        "Sample\tGenotype\tFilter\nS\tNone\tNot_assigned_to_haplotypes\n",
    );
    expect(&d.collect(), "cyp2d6.probable.diplotype", r#""*1/*2""#);
    d.put(PYP, "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\n");
    let summary = d.collect();
    expect(&summary, "cyp2d6.probable", "null");
    expect(data(&summary, "pypgx"), "low_gq", "[]");
    expect(data(&summary, "pypgx"), "not_callable", "[]");
}

#[test]
fn paraphase_counts_other_regions_and_compares_an_smn1_call() {
    let d = Sample::new("paraphase", "S");
    put_paraphase(&d, 1);
    let summary = d.collect();
    let paralogs = data(&summary, "paralogs");
    expect(section(&summary, "paralogs"), "state", r#""ok""#);
    assert_eq!(array(paralogs, "genes").len(), 4);
    for (i, gene) in ["SMN1", "SMN2", "HBA1", "HBA2"].into_iter().enumerate() {
        expect(paralogs, &format!("genes.{i}.gene"), &format!("\"{gene}\""));
    }
    expect(paralogs, "other_regions", "2");
    expect(paralogs, "other_regions_called", "1");
    expect(paralogs, "smn1.paraphase", "2");
    expect(paralogs, "smn1.parascopy.cn", "1");
    expect(paralogs, "smn1.agree", r#""disagree""#);
    put_paraphase(&d, 2);
    expect(data(&d.collect(), "paralogs"), "smn1.agree", r#""agree""#);
    fs::remove_file(d.path(SMN)).unwrap();
    expect(data(&d.collect(), "paralogs"), "smn1.parascopy", "null");
}

#[test]
fn repeat_sizing_uses_flanking_spanning_and_in_repeat_evidence() {
    let d = Sample::new("repeat-sizing", "S");
    put_repeat_sizing(&d);
    let summary = d.collect();
    let expansions = data(&summary, "expansions");
    let Value::Array(loci) = at(expansions, "loci") else {
        panic!("missing repeat loci")
    };
    assert_eq!(loci.len(), 10);
    assert_eq!(
        loci[..5]
            .iter()
            .map(|item| string(item, "locus"))
            .collect::<Vec<_>>(),
        EH_LOCI
    );
    let locus = |name: &str| {
        loci.iter()
            .find(|row| matches!(at(row, "locus"), Value::Str(s) if s == name))
            .unwrap()
    };
    expect(locus("HTT"), "normal_max", r#""26""#);
    expect(locus("HTT"), "zone", r#""27-35""#);
    assert!(string(locus("HTT"), "zone_label").contains("may expand in offspring"));
    assert!(string(locus("HTT"), "zone_citation").contains("NBK1305"));
    expect(locus("FMR1"), "zone", r#""""#);
    expect(locus("ATXN1"), "zone", r#""""#);
    expect(locus("DMPK"), "zone", r#""35-49""#);
    expect(locus("LOCUS_FL"), "sizing", r#""flanking_only""#);
    expect(locus("LOCUS_SP"), "sizing", r#""sized""#);
    expect(
        locus("LOCUS_LOW"),
        "allele_sizing",
        r#""sized/low_confidence""#,
    );
    expect(locus("LOCUS_LOW"), "sizing", r#""low_confidence""#);
    expect(locus("LOCUS_NC"), "sizing", r#""no_call""#);
    expect(locus("RFC1"), "sizing", r#""in_repeat""#);
    expect(
        expansions,
        "unsized",
        r#"["LOCUS_FL","LOCUS_LOW","LOCUS_NC","RFC1"]"#,
    );
    let Value::Array(flagged) = at(expansions, "flagged") else {
        panic!("missing flags")
    };
    let rfc1 = flagged
        .iter()
        .find(|row| matches!(at(row, "locus"), Value::Str(s) if s == "RFC1"))
        .unwrap();
    expect(rfc1, "rfc1_motif.verdict", r#""benign_motif""#);

    d.old(RFC1);
    let summary = d.collect();
    expect(
        row(data(&summary, "expansions"), "flagged", "locus", "RFC1"),
        "rfc1_motif.verdict",
        r#""not_checked""#,
    );
}

#[test]
fn annotsv_marks_exon_edges_checks_support_and_keeps_downgrades_apart() {
    const RAW: &str = "annotsv/S_sv_class4_5.tsv";
    const CHECKED: &str = "annotsv/S_sv_class4_5_checked.tsv";
    const CONSENSUS_VCF: &str = "sv_merged/S_sv_consensus.vcf.gz";
    let d = Sample::new("annotsv-support", "S");
    let raw = "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tID\tINFO\tAnnotation_mode\tGene_name\tB_loss_source\tB_loss_AFmax\tACMG_class\tAnnotSV_ranking_score\n\
        20_1000201_1009800_DEL_1\t20\t1000200\t1009800\t-9600\tDEL\tSV1\tEND=1009800;SVTYPE=DEL;SVCALLER=manta\tfull\tFIXG1\tgnomAD-SV_v4.1_DEL_fixture\t0.05\t4\t0.92\n\
        20_10101000_10102999_DEL_1\t20\t10100999\t10102999\t-2000\tDEL\tSV2\tEND=10102999;SVTYPE=DEL;SVCALLER=manta\tfull\tG1;G2;G3;G4;G5\t\t\t5\t1.20\n\
        20_20000001_20450000_DEL_1\t20\t20000000\t20450000\t-450000\tDEL\tSV3\tEND=20450000;SVTYPE=DEL;SVCALLER=manta\tfull\tG9\t\t\t4\t0.90\n";
    d.put(RAW, raw);
    d.put("genes.bed", "20\t1000000\t1010000\t+\tFIXG1\tNM_900001\t1000100\t1009900\t1000000,1006000,1009800,\t1000200,1006100,1010000,\n");
    d.put("rules.tsv", "#chrom\tstart\tend\trule\taction\tname\tsource\nchr20\t19990000\t20500000\tig_tr_locus\tclass_3\tIGFIX\tfixture\n");
    let check = |use_rules: bool| {
        let mut args = vec![
            "--class45".into(),
            d.path(RAW).display().to_string(),
            "--genes".into(),
            d.path("genes.bed").display().to_string(),
            "--out".into(),
            d.path(CHECKED).display().to_string(),
        ];
        if use_rules {
            args.extend(["--rules".into(), d.path("rules.tsv").display().to_string()]);
        }
        assert_eq!(og_steps::sv::postfilter::cli(&args), 0);
    };

    check(false);
    let consensus_header = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\n";
    d.gzip(CONSENSUS_VCF, &format!("{consensus_header}chr20\t10101000\tx\tN\t<DEL>\t.\tPASS\tSUPP=2;SUPP_VEC=11;SVTYPE=DEL;END=10103000\tGT:ID\t0/1:SV2\t0/1:DELLY2\n"));
    let summary = d.collect();
    let sv = data(&summary, "annotsv");
    expect(section(&summary, "annotsv"), "state", r#""ok""#);
    expect(sv, "count", "3");
    expect(sv, "downgraded", "[]");
    expect(sv, "rows.0.exon_edge", r#""yes""#);
    expect(sv, "rows.0.population_af", "0.05");
    expect(sv, "rows.0.common", "true");
    expect(sv, "rows.1.support", "2");

    check(true);
    let summary = d.collect();
    let sv = data(&summary, "annotsv");
    expect(sv, "count", "2");
    expect(sv, "downgraded.0.sv", r#""20:20000000-20450000""#);
    expect(sv, "downgraded.0.annotsv_class", r#""4""#);
    expect(sv, "downgraded.0.acmg_class", r#""3""#);
    expect(sv, "downgraded.0.rule", r#""ig_tr_locus""#);
    assert_eq!(array(sv, "downgraded").len(), 1);
    assert!(string(sv, "downgraded.0.rule_detail").contains("IGFIX"));
    assert!(
        array(sv, "rows")
            .iter()
            .all(|row| string(row, "sv") != "20:20000000-20450000")
    );

    check(false);
    d.gzip(CONSENSUS_VCF, consensus_header);
    let summary = d.collect();
    let sv = data(&summary, "annotsv");
    assert_eq!(
        array(sv, "rows")
            .iter()
            .map(|row| at(row, "support"))
            .collect::<Vec<_>>(),
        [&Value::Int(1); 3]
    );

    fs::remove_file(d.path(CONSENSUS_VCF)).unwrap();
    let summary = d.collect();
    assert!(
        array(data(&summary, "annotsv"), "rows")
            .iter()
            .all(|row| at(row, "support") == &Value::Null)
    );
    d.put(CONSENSUS_VCF, b"\x1f\x8b\x08\x00not gzip\n");
    let summary = d.collect();
    expect(data(&summary, "annotsv"), "count", "3");
    assert!(
        matches!(at(data(&summary, "annotsv"), "consensus_note"), Value::Str(s) if s.contains("could not be read"))
    );
    d.put(CHECKED, raw);
    expect(section(&d.collect(), "annotsv"), "state", r#""unreadable""#);
    fs::remove_file(d.path(CHECKED)).unwrap();
    let summary = d.collect();
    expect(section(&summary, "annotsv"), "state", r#""missing""#);
    assert!(
        array(&summary, "not_run")
            .iter()
            .any(|item| matches!(item, Value::Str(s) if s.contains("AnnotSV")))
    );
}

fn array<'a>(value: &'a Value, path: &str) -> &'a [Value] {
    let Value::Array(items) = at(value, path) else {
        panic!("{path} is not an array")
    };
    items
}

fn string<'a>(value: &'a Value, path: &str) -> &'a str {
    let Value::Str(text) = at(value, path) else {
        panic!("{path} is not a string")
    };
    text
}

fn row<'a>(value: &'a Value, path: &str, key: &str, name: &str) -> &'a Value {
    array(value, path)
        .iter()
        .find(|entry| string(entry, key) == name)
        .unwrap_or_else(|| panic!("missing {name} in {path}"))
}

#[test]
fn stranger_flags_full_mutations_and_ignores_its_outdated_vcf() {
    let d = Sample::new("stranger-status", "S");
    put_stranger(&d);
    let summary = d.collect();
    expect(
        section(&summary, "expansions"),
        "source",
        r#""expansion_hunter/S_eh_stranger.vcf""#,
    );
    let expansions = data(&summary, "expansions");
    expect(expansions, "stranger", "true");
    expect(expansions, "no_status", "1");
    expect(
        expansions,
        "flagged",
        r#"[{"locus":"LOCUS_B","repeat_count":"20/75","status":"full_mutation","sizing":"unknown"}]"#,
    );
    assert_eq!(
        array(expansions, "key_loci")
            .iter()
            .map(|item| string(item, "locus"))
            .collect::<Vec<_>>(),
        EH_LOCI
    );

    d.old(STRANGER);
    d.put(STR, fs::read(d.path(STR)).unwrap());
    let summary = d.collect();
    expect(
        section(&summary, "expansions"),
        "source",
        r#""expansion_hunter/S_eh.vcf""#,
    );
    let expansions = data(&summary, "expansions");
    expect(expansions, "stranger", "false");
    expect(expansions, "stranger_outdated", "true");
}

fn repeat_reading<'a>(summary: &'a Value, locus: &str) -> &'a Value {
    row(data(summary, "expansions"), "readings", "locus", locus)
}

fn ugt1a1_call(summary: &Value, genotype: &str, agreement: &str) {
    let ugt = repeat_reading(summary, "UGT1A1");
    assert_eq!(string(ugt, "genotype"), genotype);
    assert_eq!(string(ugt, "agreement"), agreement);
}

#[test]
fn extra_repeats_keep_the_two_cftr_tracts_and_ugt1a1_agreement() {
    let d = Sample::new("extra-repeats-exa", "S");
    put_extra_repeats(&d, extra_repeats_exa());
    let summary = d.collect();
    let expansions = data(&summary, "expansions");
    assert_eq!(
        array(expansions, "loci")
            .iter()
            .map(|locus| string(locus, "locus"))
            .filter(|locus| locus.starts_with("CFTR"))
            .collect::<Vec<_>>(),
        ["CFTR_T", "CFTR_TG"]
    );
    expect(
        expansions,
        "not_graded",
        r#"{"CFTR":"the CFTR tract","UGT1A1":"the UGT1A1 repeat"}"#,
    );
    let cftr = repeat_reading(&summary, "CFTR");
    expect(cftr, "genotype", r#""TG11/TG12 5T/7T""#);
    let pairings = array(cftr, "pairings");
    assert_eq!(pairings.len(), 2);
    assert!(pairings.contains(&Value::Str("TG12-5T/TG11-7T".into())));
    expect(cftr, "five_t", "true");
    expect(cftr, "phase", r#""not_needed""#);
    assert!(string(cftr, "reading").starts_with("5T: 1 copy with the 5T poly-T tract"));
    assert!(string(cftr, "reading").contains("No R117H copy"));
    ugt1a1_call(&summary, "*1/*28", "agrees");
    assert!(string(repeat_reading(&summary, "UGT1A1"), "reading").contains("TA6/TA7"));
    let slc6a4 = repeat_reading(&summary, "SLC6A4");
    expect(slc6a4, "status", r#""not_tested""#);
    assert!(string(slc6a4, "source").contains("37032427"));
}

#[test]
fn extra_repeats_without_5t_or_a_pharmcat_call() {
    let d = Sample::new("extra-repeats-exb", "S");
    put_extra_repeats(&d, extra_repeats_exb());
    let summary = d.collect();
    let cftr = repeat_reading(&summary, "CFTR");
    expect(cftr, "genotype", r#""TG11-7T/TG11-7T""#);
    expect(cftr, "five_t", "false");
    assert!(string(cftr, "reading").starts_with("No 5T copy"));
    ugt1a1_call(&summary, "*1/*1", "not_compared");
}

#[test]
fn extra_repeats_keep_unknown_phase_and_pharmcat_disagreement() {
    let d = Sample::new("extra-repeats-exc", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            tg: "12/12",
            t: "5/7",
            ugt: "8/8",
            r117h: Some("G/A"),
            pharmcat: Some("*1/*1"),
            ..ExtraRepeats::default()
        },
    );
    let summary = d.collect();
    let cftr = repeat_reading(&summary, "CFTR");
    expect(cftr, "genotype", r#""TG12-5T/TG12-7T""#);
    expect(cftr, "phase", r#""unknown""#);
    assert!(string(cftr, "reading").contains("phase unknown"));
    assert!(string(cftr, "reading").contains("TG12 on the 5T copy"));
    ugt1a1_call(&summary, "*28/*28", "differs");
}

#[test]
fn extra_repeats_put_r117h_in_cis_with_two_5t_copies_and_name_ugt1a1_stars() {
    let d = Sample::new("extra-repeats-exd", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            tg: "11/12",
            t: "5/5",
            ugt: "6/9",
            r117h: Some("G/A"),
            ..ExtraRepeats::default()
        },
    );
    let summary = d.collect();
    expect(repeat_reading(&summary, "CFTR"), "phase", r#""in_cis""#);
    assert_eq!(
        string(repeat_reading(&summary, "UGT1A1"), "genotype"),
        "*36/*37"
    );
}

#[test]
fn extra_repeats_without_opt_in_explain_how_to_add_the_loci() {
    let d = Sample::new("extra-repeats-exe", "S");
    put_extra_repeats(&d, ExtraRepeats::default());
    let summary = d.collect();
    for locus in ["CFTR", "UGT1A1"] {
        let reading = repeat_reading(&summary, locus);
        expect(reading, "status", r#""not_sized""#);
        assert!(
            string(reading, "reading")
                .contains("STR_EXTRA_LOCI=yes ./scripts/setup.sh --strchive-data <genome_dir>")
        );
    }
    expect(
        repeat_reading(&summary, "SLC6A4"),
        "status",
        r#""not_tested""#,
    );
}

fn combined_ugt1a1_agrees(case: &str, pharmcat: &str) {
    let d = Sample::new(case, "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            ugt: "7/8",
            pharmcat: Some(pharmcat),
            ..ExtraRepeats::default()
        },
    );
    ugt1a1_call(&d.collect(), "*1/*28", "agrees");
}

#[test]
fn extra_repeats_compare_a_combined_pharmcat_allele() {
    combined_ugt1a1_agrees("extra-repeats-exf", "*80+*28/*1");
}

#[test]
fn extra_repeats_compare_a_bracketed_pharmcat_allele() {
    combined_ugt1a1_agrees("extra-repeats-exg", "[*80 + *28]/*1");
}

#[test]
fn extra_repeats_describe_one_r117h_copy_without_5t() {
    let d = Sample::new("extra-repeats-exh", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            t: "7/9",
            r117h: Some("G/A"),
            ..extra_repeats_exb()
        },
    );
    let summary = d.collect();
    let cftr = repeat_reading(&summary, "CFTR");
    expect(cftr, "five_t", "false");
    expect(cftr, "r117h", r#""G/A""#);
    let reading = string(cftr, "reading");
    assert!(reading.contains("R117H (G/A) is read"));
    assert!(reading.contains("7T or 9T copy"));
    assert!(!reading.contains("does not change how a CFTR finding is read"));
}

#[test]
fn extra_repeats_without_a_traits_table_say_r117h_was_not_read() {
    let d = Sample::new("extra-repeats-exi", "S");
    put_extra_repeats(&d, extra_repeats_exb());
    assert!(string(repeat_reading(&d.collect(), "CFTR"), "reading").contains("R117H was not read"));
}

#[test]
fn extra_repeats_do_not_assign_function_to_an_unnamed_ta_length() {
    let d = Sample::new("extra-repeats-exj", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            ugt: "7/11",
            pharmcat: Some("*1/*1"),
            ..ExtraRepeats::default()
        },
    );
    let summary = d.collect();
    ugt1a1_call(&summary, "*1/TA10", "not_compared");
    let reading = string(repeat_reading(&summary, "UGT1A1"), "reading");
    assert!(reading.contains("not compared with PharmCAT"));
    assert!(!reading.contains("decreased-function"));
}

#[test]
fn extra_repeats_put_two_r117h_copies_on_both_poly_t_lengths() {
    let d = Sample::new("extra-repeats-exk", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            t: "7/9",
            r117h: Some("A/A"),
            ..extra_repeats_exb()
        },
    );
    let summary = d.collect();
    let cftr = repeat_reading(&summary, "CFTR");
    expect(cftr, "five_t", "false");
    let reading = string(cftr, "reading");
    assert!(reading.contains("R117H (A/A) is read: it sits on both copies (7T and 9T)"));
    assert!(!reading.contains("it sits on a 7T or 9T copy"));
}

#[test]
fn extra_repeats_put_two_r117h_copies_on_both_7t_copies() {
    let d = Sample::new("extra-repeats-exl", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            r117h: Some("A/A"),
            ..extra_repeats_exb()
        },
    );
    assert!(
        string(repeat_reading(&d.collect(), "CFTR"), "reading")
            .contains("R117H (A/A) is read: it sits on both 7T copies")
    );
}

#[test]
fn extra_repeats_explain_why_an_unreadable_traits_table_cannot_phase_r117h() {
    let d = Sample::new("extra-repeats-bad-traits", "S");
    put_extra_repeats(
        &d,
        ExtraRepeats {
            tg: "12/12",
            t: "5/7",
            ugt: "7/8",
            ..ExtraRepeats::default()
        },
    );
    d.put(TRAITS, "rsid\tgenotype\tstatus\nrs78655421\tG/A\tcalled\n");
    let summary = d.collect();
    let cftr = repeat_reading(&summary, "CFTR");
    expect(cftr, "phase", r#""not_read""#);
    assert!(string(cftr, "reading").contains("traits table is unreadable"));
}

#[test]
fn extra_repeats_identify_traits_and_pharmcat_tables_from_earlier_runs() {
    let d = Sample::new("extra-repeats-stale-side-tables", "S");
    let status = format!(
        "meta\tstarted_epoch\t{STARTED}\nstep\t09\tok\nstep\t27\tskipped\nstep\t42\tskipped\n"
    );
    put_extra_repeats(
        &d,
        ExtraRepeats {
            tg: "12/12",
            t: "5/7",
            ugt: "7/8",
            r117h: Some("G/A"),
            pharmcat: Some("*1/*28"),
            side_epoch: Some(STARTED - 60),
            run_status: Some(&status),
            ..ExtraRepeats::default()
        },
    );
    let summary = d.collect();
    assert!(
        string(repeat_reading(&summary, "CFTR"), "reading")
            .contains("R117H from step 42's traits table of an earlier run")
    );
    ugt1a1_call(&summary, "*1/*28", "agrees");
    assert!(
        string(repeat_reading(&summary, "UGT1A1"), "reading")
            .contains("step 27's table is from an earlier run")
    );
}

#[test]
fn mito_polyc_and_allele_specific_annotation_keep_distinct_readings() {
    let d = Sample::new("mito-annotation", "S");
    put_mito_calls(&d);
    let summary = d.collect();
    let mito = data(&summary, "mito");
    expect(mito, "pass", "4");
    expect(mito, "heteroplasmic", "1");
    assert_eq!(
        array(mito, "polyc_heteroplasmic")
            .iter()
            .map(|site| string(site, "variant"))
            .collect::<Vec<_>>(),
        ["chrM:302 AC>ACC,A", "chrM:16189 T>C"]
    );
    expect(
        mito,
        "heteroplasmic_sites",
        r#"[{"pos":3243,"ref":"A","alt":"G","af":0.3,"origin_edge":false,"low_af":false}]"#,
    );

    put_mito_databases(&d);
    annotate_mito(&d, true);
    let table = fs::read_to_string(d.path(MITO_ANNOTATION)).unwrap();
    let mut lines = table.lines().filter(|line| !line.starts_with("# "));
    let headers: Vec<_> = lines.next().unwrap().split('\t').collect();
    let rows: Vec<std::collections::BTreeMap<_, _>> = lines
        .map(|line| headers.iter().copied().zip(line.split('\t')).collect())
        .collect();
    assert_eq!(rows.len(), 6);
    let allele = |pos, reference, alt| {
        rows.iter()
            .find(|row| (row["pos"], row["ref"], row["alt"]) == (pos, reference, alt))
            .unwrap()
    };
    let insertion = allele("302", "A", "AC");
    assert_eq!(
        (insertion["gnomad_af_het"], insertion["af"]),
        ("0.2", "0.4")
    );
    assert_eq!(allele("302", "AC", "A")["af"], "0.1");
    let hit = allele("3243", "A", "G");
    assert_eq!(
        (
            hit["clinvar_significance"],
            hit["clinvar_review"],
            hit["clinvar_condition"],
            hit["gnomad_af_hom"],
            hit["gnomad_af_het"]
        ),
        (
            "Pathogenic",
            "reviewed_by_expert_panel",
            "MELAS",
            "0.00002",
            "0.0003"
        )
    );
    assert_eq!(allele("16189", "T", "C")["clinvar_significance"], ".");

    let summary = d.collect();
    let annotation = at(data(&summary, "mito"), "annotation");
    expect(annotation, "clinvar_plp", "1");
    assert!(
        array(annotation, "variants")
            .iter()
            .all(|item| !matches!(string(item, "pos"), "73" | "9000"))
    );

    annotate_mito(&d, false);
    expect(data(&d.collect(), "mito"), "annotation.clinvar_plp", "0");
}

#[test]
fn parascopy_catalog_survives_missing_smn_and_unreadable_loci() {
    let d = Sample::new("parascopy-catalog", "S");
    put_parascopy_catalog(&d);
    let summary = d.collect();
    expect(section(&summary, "parascopy"), "state", r#""ok""#);
    let p = data(&summary, "parascopy");
    expect(p, "hba.alpha_globin_genes", r#""4""#);
    assert!(!hba_not_assessed(&summary));
    expect(p, "smn1.cn", "null");
    expect(
        p,
        "reading",
        r#""no SMN1 table (step 35 wrote no SMN copy-number table): no carrier reading""#,
    );
    assert_eq!(
        array(p, "loci")
            .iter()
            .map(|locus| string(locus, "table"))
            .collect::<Vec<_>>(),
        ["gba1", "c4", "amy1", "mine"]
    );
    let locus = |name| row(p, "loci", "table", name);
    let gba = locus("gba1");
    expect(gba, "cn", "4");
    expect(gba, "status", r#""called""#);
    expect(
        gba,
        "text",
        r#""4 (2 region(s), quality 38 or more, filter PASS)""#,
    );
    expect(gba, "genes", r#""GBA1,GBAP1""#);
    assert!(string(gba, "reading").contains("Gaucher"));
    let c4 = locus("c4");
    expect(c4, "cn", "null");
    expect(c4, "status", r#""called""#);
    expect(c4, "reading", r#"".""#);
    assert!(string(c4, "text").contains("chr6:31980532-31990000 5 (quality 40, filter PASS)"));
    assert!(string(c4, "text").contains("4 (quality 12, filter LowQual)"));
    expect(locus("amy1"), "status", r#""not assessed""#);
    expect(
        locus("amy1"),
        "text",
        r#""not assessed: no read over chr1:103745407-103761087""#,
    );
    let mine = locus("mine");
    expect(mine, "cn", "2");
    expect(mine, "genes", r#""mine""#);
    expect(mine, "reading", r#"".""#);

    d.put(C4, "chrom\tstart\tend\tagCN\nchr6\t1\t2\t4\n");
    let summary = d.collect();
    let p = data(&summary, "parascopy");
    expect(row(p, "loci", "table", "c4"), "status", r#""unreadable""#);
    expect(p, "hba.alpha_globin_genes", r#""4""#);

    d.put(
        "paralogs/S_notch2_copy_number.tsv",
        format!("{SMN_HEADER}chr1\t120000000\t120100000\tNOTCH2\tPASS\t<5\t30\t.\t.\t.\t*\n"),
    );
    let summary = d.collect();
    let bound = row(data(&summary, "parascopy"), "loci", "table", "notch2");
    expect(bound, "status", r#""not resolved""#);
    expect(bound, "cn", "null");
    assert!(
        string(bound, "text")
            .ends_with(": no region with a whole copy number, filter PASS and quality 20 or more")
    );
}

fn hba_not_assessed(summary: &Value) -> bool {
    array(summary, "not_assessed")
        .iter()
        .any(|item| matches!(item, Value::Str(text) if text.starts_with("HBA1/HBA2 copy number")))
}

#[test]
fn hba_is_assessed_by_depth_or_by_two_called_paraphase_genes() {
    let d = Sample::new("hba-assessment", "S");
    d.put(SMN, format!("{SMN_HEADER}{SMN_ELSEWHERE}{SMN_EXON_7}"));
    let summary = d.collect();
    assert!(hba_not_assessed(&summary));
    expect(
        data(&summary, "parascopy"),
        "smn1_marker.status",
        r#""not assessed""#,
    );
    d.put(MARKER, "metric\tvalue\nstatus\tmaybe\n");
    expect(
        data(&d.collect(), "parascopy"),
        "smn1_marker.status",
        r#""not assessed""#,
    );

    d.put(HBA, HBA_DEPTH_3);
    let summary = d.collect();
    expect(
        data(&summary, "parascopy"),
        "hba.alpha_globin_genes",
        r#""3""#,
    );
    assert!(!hba_not_assessed(&summary));
    d.put(HBA, HBA_LOW);
    assert!(hba_not_assessed(&d.collect()));

    for path in [SMN, MARKER, HBA] {
        fs::remove_file(d.path(path)).unwrap();
    }
    let summary = d.collect();
    for step in [
        "SMN1/SMN2 copy number (Parascopy, short reads)",
        "Paralog genes (Paraphase, PacBio HiFi)",
    ] {
        assert!(array(&summary, "not_run").contains(&Value::Str(step.into())));
    }
    for (rows, not_assessed) in [
        (format!("{HBA1}{HBA2}"), false),
        ("smn1\tSMN1\t2\tcalled\t.\t4\t90\t31\n".into(), true),
        (HBA1.into(), true),
        (format!("{HBA1}hba\tHBA2\t.\tno call\t.\t.\t90\t31\n"), true),
    ] {
        d.put(PARAPHASE_CN, format!("{PARAPHASE_HEADER}{rows}"));
        assert_eq!(hba_not_assessed(&d.collect()), not_assessed, "{rows}");
    }
}

#[test]
fn cn_tables_read_dosage_only_and_reject_cut_short_rows_and_headers() {
    let only = Sample::new("dosage-only", "S");
    put_dosage_only(&only);
    let summary = only.collect();
    let dosage = array(data(&summary, "parascopy"), "dosage");
    assert_eq!(
        dosage
            .iter()
            .map(|row| string(row, "copies"))
            .collect::<Vec<_>>(),
        ["3"]
    );

    let d = Sample::new("cn-tables", "S");
    d.put(SMN, format!("{SMN_HEADER}{SMN_EXON_7}"));
    put_cn_tables(&d);
    let summary = d.collect();
    let p = data(&summary, "parascopy");
    expect(p, "cyp21a2.finding", r#""deletion""#);
    assert_eq!(
        array(p, "dosage")
            .iter()
            .map(|row| (
                string(row, "locus"),
                string(row, "status"),
                string(row, "copies")
            ))
            .collect::<Vec<_>>(),
        [
            ("HP", "called", "3"),
            ("LCE3BC", "called", "1"),
            ("CCL3L1", "not assessed", ".")
        ],
    );

    for (path, contents) in [
        (DOSAGE, format!("{DOSAGE_HEADER}HP\tcalled\n")),
        (
            INTERVALS,
            format!("{IV_COLS}{INTERVAL_UNCHANGED}AZFc\tAZFc\n"),
        ),
        (DOSAGE, "locus\tstatus\n".into()),
        (INTERVALS, "locus\tname\tstatus\n".into()),
        (SMN, "chrom\tstart\tend\tagCN\n".into()),
    ] {
        d.put(path, &contents);
        expect(
            section(&d.collect(), "parascopy"),
            "state",
            r#""unreadable""#,
        );
        fs::remove_file(d.path(path)).unwrap();
    }

    d.put(SMN, format!("{SMN_HEADER}{SMN_EXON_7}"));
    for (header, status) in [
        ("chrom\tstart\tend\tlocus\n", r#""unreadable""#),
        (SMN_HEADER, r#""not assessed""#),
    ] {
        d.put(GBA, header);
        let summary = d.collect();
        expect(
            row(data(&summary, "parascopy"), "loci", "table", "gba1"),
            "status",
            status,
        );
    }
}

#[test]
fn hla_keeps_single_allele_quality_and_only_withholds_genes_pharmcat_uses() {
    let d = Sample::new("hla-allele-quality", "S");
    put_hla_allele_quality(&d);
    let summary = d.collect();
    let hla = data(&summary, "hla");
    let locus = |name| row(hla, "loci", "gene", name);
    expect(locus("HLA-A"), "alleles", r#"["A*01:01:01"]"#);
    expect(locus("HLA-A"), "quality", r#"["60"]"#);
    expect(locus("HLA-A"), "low_confidence", "false");
    expect(locus("HLA-B"), "quality", r#"["0","40"]"#);
    expect(locus("HLA-B"), "low_confidence", "true");
    expect(locus("HLA-B"), "withheld_from_pharmcat", "true");
    expect(locus("HLA-C"), "low_confidence", "true");
    expect(locus("HLA-C"), "withheld_from_pharmcat", "false");
    let Value::Object(fields) = hla else {
        panic!("missing HLA data")
    };
    for key in ["optitype", "kir"] {
        assert!(matches!(fields.get(key), None | Some(Value::Null)), "{key}");
    }
}

#[test]
fn roh_upd_and_relationship_depend_on_the_distribution_of_long_runs() {
    let catalog =
        fs::read_to_string(std::path::Path::new(REPO).join("assets/imprinted_regions.tsv"))
            .unwrap();
    let rows: Vec<Vec<_>> = catalog
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        rows.iter().map(|row| row[0]).collect::<Vec<_>>(),
        ["chr6", "chr7", "chr11", "chr14", "chr15", "chr20"]
    );
    assert!(
        rows.iter()
            .all(|row| row[6].split(';').all(|pmid| pmid.parse::<u64>().is_ok()))
    );

    let d = Sample::new("roh-upd", "S");
    put_roh_with_sites(&d, UPD_RUNS);
    let summary = d.collect();
    let roh = data(&summary, "roh");
    expect(roh, "upd_flag.chrom", r#""chr15""#);
    expect(roh, "upd_flag.mb", "81.5");
    expect(roh, "relationship", "null");
    expect(roh, "upd_flag.imprinted.0.maternal_upd", r#""-""#);
    expect(roh, "upd_flag.imprinted.0.paternal_upd", r#""-""#);
    check_f_roh(roh, ROH_UPD_LONG_BP);

    for (runs, chr, region) in [
        (UPD_RUNS, "chr15", "15q11-q13 (SNRPN, UBE3A)"),
        (UPD_DISTAL, "chr15", "15q11-q13 (SNRPN, UBE3A)"),
        (
            "RG\tS\tchr7\t70000001\t110000000\t40000000\t900\t50\n",
            "chr7",
            "7p12 (GRB10), 7q32 (MEST)",
        ),
    ] {
        put_roh_with_sites(&d, runs);
        let summary = d.collect();
        let roh = data(&summary, "roh");
        assert_eq!(string(roh, "upd_flag.chrom"), chr);
        let imprinted = array(roh, "upd_flag.imprinted");
        assert_eq!(imprinted.len(), 1);
        assert_eq!(string(&imprinted[0], "region"), region);
        assert_eq!(string(&imprinted[0], "maternal_upd"), "-");
    }

    put_roh_with_sites(&d, UPD_OVER_HALF);
    let summary = d.collect();
    expect(data(&summary, "roh"), "upd_flag.chrom", r#""chr2""#);

    put_roh_spread(&d);
    let summary = d.collect();
    let roh = data(&summary, "roh");
    expect(roh, "upd_flag", "null");
    check_f_roh(roh, 12 * SPREAD_LENGTH);
    expect(
        roh,
        "relationship",
        r#""first cousins once removed (F 1/32)""#,
    );

    put_roh_with_sites(&d, "RG\tS\tchr2\t10000001\t13000000\t3000000\t900\t50\n");
    let summary = d.collect();
    let roh = data(&summary, "roh");
    check_f_roh(roh, 0);
    expect(roh, "relationship", r#""no sign of related parents""#);

    put_roh_with_sites(&d, "RG\tS\tchr2\t10000001\t15000000\t5000000\t900\t50\n");
    let summary = d.collect();
    let roh = data(&summary, "roh");
    expect(roh, "long_mb", "5.0");
    assert_eq!(array(roh, "autosomal_over_5mb").len(), 1);
}

fn check_f_roh(roh: &Value, bases: u64) {
    let Value::Float(value) = at(roh, "f_roh") else {
        panic!("F_ROH is not numeric")
    };
    assert!((value - bases as f64 / ROH_DENOM).abs() < 0.001, "{value}");
}

#[test]
fn cnvpytor_counts_significant_calls_and_filters_by_sex_and_quality() {
    let d = Sample::new("cnvpytor-filters", "S");
    put_cnv_calls(&d);
    for (sex, kept) in [(Some("female"), "1"), (Some("male"), "2"), (None, "2")] {
        let summary = d.collect_with_sex(sex);
        let cnv = data(&summary, "cnvpytor");
        expect(cnv, "significant", "5");
        expect(cnv, "significant_filtered", kept);
        if sex == Some("female") {
            expect(
                cnv,
                "dropped",
                r#"{"q0":1,"gap":1,"band":1,"chry_female":1}"#,
            );
        }
    }
}

#[test]
fn karyotypes_use_plausible_copy_numbers_and_enough_x_sites() {
    let d = Sample::new("karyotype-copy-numbers", "S");
    for (x, y, name) in KARYOTYPES {
        put_karyotype(&d, x, y);
        let summary = d.collect();
        let want = name.map_or(Value::Null, |name| Value::Str(name.into()));
        for source in ["sex_check", "sample_qc"] {
            assert_eq!(
                at(data(&summary, source), "karyotype"),
                &want,
                "{source}: {x}, {y}"
            );
        }
    }
    for sites in ["10", "", "inf"] {
        put_somalier_ratios(&d, sites, 1.02, 0.49);
        expect(data(&d.collect(), "sample_qc"), "karyotype", "null");
    }
}

#[test]
fn acmg_secondary_findings_list_only_panel_genes_and_high_impact() {
    let d = Sample::new("acmg-secondary-findings", "S");
    put_acmg_secondary_findings(&d);
    let summary = d.collect();
    let sf = at(data(&summary, "clinical"), "acmg_sf");
    expect(sf, "version", r#""ACMG SF v3.3""#);
    expect(sf, "genes_on_list", "84");
    let clinvar = array(sf, "clinvar_hits");
    assert_eq!(clinvar.len(), 1);
    expect(&clinvar[0], "gene", r#""BRCA2""#);
    let high_impact = array(sf, "high_impact");
    assert_eq!(high_impact.len(), 1);
    expect(&high_impact[0], "gene", r#""TTN""#);
    expect(&high_impact[0], "consequence", r#""stop_gained""#);
}

#[test]
fn prs_format_copies_the_legacy_score_evaluations_and_applies_to() {
    use og_steps::prs::{self, Args};
    use std::io::Read;

    let d = Sample::new("prs-format-evaluations", "S");
    d.gzip("scores/PGS000018.txt.gz", "#pgs_id=PGS000018\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t1100\tA\tG\t0.5\n");
    d.put(
        "plain.tsv",
        "pgs_id\ttrait_reported\nPGS000018\tCoronary artery disease\n",
    );
    d.put("female.tsv", "pgs_id\ttrait_reported\tapplies_to\tppm_id\tancestry\teffect\nPGS000018\tCoronary artery disease\tfemale\tPPM000591\tEUR\tHR 1.53\n");
    d.put(
        "results/sample/match/sample_summary.csv",
        "dataset,accession,match_status,count\nsample,PGS000018,matched,1\n",
    );
    d.gzip("results/sample/score/sample_pgs.txt.gz", "sampleset\tFID\tIID\tPGS\tSUM\tpercentile_MostSimilarPop\nsample\tHG\tHG\tPGS000018\t0.5\t81.2\n");
    d.gzip("results/sample/score/sample_popsimilarity.txt.gz", "sampleset\tFID\tIID\tPC1\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\nsample\tHG\tHG\t0.1\tEAS\tFalse\tFalse\n");
    let shipped = std::path::Path::new(REPO).join("assets/pgs_scores.tsv");
    for (label, applies_to, evaluations) in [
        (
            shipped,
            "",
            "PGS000018\tPPM000591;PPM000594;PPM000597\tEUR;AFR;AMR\tHR 1.53;HR 1.27;HR 1.53\n",
        ),
        (d.path("plain.tsv"), "", ""),
        (
            d.path("female.tsv"),
            "female",
            "PGS000018\tPPM000591\tEUR\tHR 1.53\n",
        ),
    ] {
        let out = d.path("formatted");
        prs::run(
            &Args {
                scores: d.path("scores"),
                labels: Some(label),
                out: out.clone(),
                alleles: d.path("alleles.tsv"),
                ..Args::default()
            },
            &mut Vec::new(),
        )
        .unwrap();
        assert_eq!(fs::read_dir(&out).unwrap().count(), 4);
        let mut score = String::new();
        flate2::read::GzDecoder::new(fs::File::open(out.join("PGS000018.txt.gz")).unwrap())
            .read_to_string(&mut score)
            .unwrap();
        assert_eq!(
            score,
            "#pgs_id=PGS000018\n#pgs_name=PGS000018\n#trait_reported=Coronary artery disease\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t1100\tA\tG\t0.5\n"
        );
        assert_eq!(
            fs::read_to_string(out.join(prs::CHRX_FILE)).unwrap(),
            "pgs_id\tchrx_rows\tchrx\nPGS000018\t0\tleft out\n"
        );
        assert_eq!(
            fs::read_to_string(out.join(prs::SCORE_NOTES_FILE)).unwrap(),
            format!(
                "pgs_id\tapplies_to\tambiguous_rows\tambiguous\nPGS000018\t{applies_to}\t0\tdropped\n"
            )
        );
        assert_eq!(
            fs::read_to_string(out.join(prs::EVALS_FILE)).unwrap(),
            format!("pgs_id\tppm_id\tancestry\teffect\n{evaluations}")
        );
        if applies_to.is_empty() && !evaluations.is_empty() {
            // The other PRS table scenarios are owned by og-summary/src/prs/tests.rs.
            // Exercise the legacy EAS case through the table verb as well as the reports.
            let mut args = vec![
                "--sample".into(),
                "S".into(),
                "--sampleset".into(),
                "sample".into(),
                "--panel".into(),
                "pgsc_1000G_v1".into(),
                "--input-kind".into(),
                "gvcf".into(),
            ];
            for (flag, path) in [
                ("--results", "results"),
                ("--scores", "formatted"),
                ("--out", PRS),
            ] {
                args.extend([flag.into(), d.path(path).display().to_string()]);
            }
            fs::create_dir_all(d.path("prs")).unwrap();
            let mut err = Vec::new();
            assert_eq!(
                og_summary::prs::cli(&args, &mut Vec::new(), &mut err),
                0,
                "{}",
                String::from_utf8_lossy(&err)
            );
            let table = fs::read_to_string(d.path(PRS)).unwrap();
            let mut lines = table.lines();
            let header: Vec<_> = lines.next().unwrap().split('\t').collect();
            assert_eq!(&header[12..], ["Evaluated_In", "Effect", "Input"]);
            let cells: std::collections::BTreeMap<_, _> = header
                .into_iter()
                .zip(lines.next().unwrap().split('\t'))
                .collect();
            assert_eq!(
                (
                    cells["Ancestry_Group"],
                    cells["Evaluated_In"],
                    cells["Effect"]
                ),
                ("EAS", "EUR;AFR;AMR", "")
            );
        }
    }
}
