//! The collector cases from tests/test_collect_summary.py. The section map
//! lives in the PR alongside the cases already covered by owning crates.

#[path = "port_tests_summary_support/port_tests_summary_fixtures.rs"]
mod fixtures;

use fixtures::*;
use pycompat::json::{Dump, Value, dumps, loads};
use std::fs;
use std::io::Write;
use std::num::NonZeroUsize;

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
    fs::remove_file(stale.path(RUN)).unwrap();
    expect(section(&stale.collect(), "slivar"), "state", r#""ok""#);

    let failed = Sample::new("failed", "S");
    put_failed(&failed);
    let summary = failed.collect();
    expect(section(&summary, "telomere"), "state", r#""failed""#);
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
    let vcf = VCF.replacen("#CHROM", "##PGP_CALL_FLOOR=GQ=30,DP=10\n#CHROM", 1);
    let hits = format!(
        "{vcf}chr1\t100\t.\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\tGT:GQ:DP\t0/1:12:6\n\
        chr1\t40000\t.\tG\tA\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT:GQ:DP\t0/1:50:40\n\
        chr2\t20000\t.\tCA\tC\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\tGT:GQ:DP\t0/1:45:35\n\
        chr2\t40000\t.\tA\tT\t50\tPASS\tGENEINFO=GENED:44;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n"
    );
    d.put(CLINVAR, &hits);
    let path = "vep/S_vep.vcf.gz";
    let tbi = "vep/S_vep.vcf.gz.tbi";
    let info = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Existing_variation|gnomADe_AF|gnomADg_AF\">";
    let chr1 = format!(
        "{}{info}\n#CHROM{}",
        VCF.split_once("#CHROM").unwrap().0,
        VCF.split_once("#CHROM").unwrap().1
    );
    let mut bgzf = og_io::bgzf::Writer::new(Vec::new(), NonZeroUsize::MIN).unwrap();
    bgzf.write_all(chr1.as_bytes()).unwrap();
    bgzf.write_all(b"chr1\t100\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|GENEA|rs123&COSV1|0.03|0.02,T|intron_variant|MODIFIER|GENEZ|rs123|.|.\tGT\t0/1\n\
        chr1\t40000\t.\tG\tA,C\t50\tPASS\tCSQ=A|stop_gained|HIGH|GENEB|rs77&rs9|0.0002|,C|missense_variant|MODERATE|GENEB|rs555|0.9|\tGT\t0/1\n").unwrap();
    bgzf.flush().unwrap();
    bgzf.write_all(b"chr2\t20000\t.\tCA\tC\t50\tPASS\tCSQ=-|frameshift_variant|HIGH|GENEC|COSV9|0.0001|0.00005\tGT\t0/1\n").unwrap();
    let bytes = bgzf.finish().unwrap();
    let index = og_io::tabix_write::build(bytes.as_slice()).unwrap();
    d.put(path, bytes);
    index.write_tbi(&d.path(tbi), false).unwrap();

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
    for (gene, af, rsids, low, common) in [
        ("GENEA", "0.03", r#"["rs123"]"#, "true", "true"),
        ("GENEB", "0.0002", r#"["rs9","rs77"]"#, "false", "false"),
        ("GENEC", "0.0001", "[]", "false", "false"),
        ("GENED", "null", "[]", "false", "false"),
    ] {
        let hit = found
            .iter()
            .find(|hit| at(hit, "gene") == &Value::Str(gene.into()))
            .unwrap();
        expect(hit, "af", af);
        expect(hit, "rsids", rsids);
        expect(hit, "low_confidence", low);
        expect(hit, "common", common);
    }

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

    d.put(CLINVAR, hits.replace("GQ=30,DP=10", "GQ=0,DP=0"));
    expect(data(&d.collect(), "clinvar"), "count", "4");
    expect(data(&d.collect(), "clinvar"), "low_confidence", "0");
    fs::remove_file(d.path(tbi)).unwrap();
    let no_index = d.collect();
    let cv = data(&no_index, "clinvar");
    assert!(matches!(at(cv, "frequency.note"), Value::Str(s) if s.contains("no .tbi")));
    fs::remove_file(d.path(path)).unwrap();
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
    expect(&summary, "databases.cpsr_clinvar_release", r#""20250309""#);

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
    expect(couple, "partner", r#""A""#);
    expect(couple, "shared.0.gene", r#""NAGA""#);
    expect(couple, "shared.0.own", r#""b-snv""#);
    expect(couple, "shared.0.partner_findings", r#""a-snv""#);
    expect(couple, "x_linked.0.gene", r#""FMR1""#);
    expect(couple, "one_partner", "1");
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
    d.old(CONSENSUS);
    expect(&d.collect(), "cyp2d6.consensus", "null");
}

#[test]
fn prs_keeps_the_group_where_the_effect_was_evaluated() {
    for (group, effect) in [("EUR", "HR 1.53 (PPM000591)"), ("EAS", "")] {
        let d = Sample::new(group, "S");
        put_prs_evaluation(&d, group, effect);
        let summary = d.collect();
        let row = at(data(&summary, "prs"), "scores.0");
        expect(row, "evaluated_in", r#"["EUR"]"#);
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
    expect(data(&no_catalog, "expansions"), "no_status", "2");
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
    expect(section(&summary, "parascopy"), "step", r#""35""#);
    expect(p, "smn1.cn", "2");
    expect(p, "smn2.cn", "2");
    expect(p, "agcn", r#""4""#);
    expect(p, "agcn_qual", r#""60""#);
    expect(p, "hba.alpha_globin_genes", r#""4""#);
    expect(p, "smn1_marker.status", r#""present""#);

    d.put(SMN, format!("{SMN_HEADER}chr5\t70075495\t70077900\tSMN1\tPASS\t3\t50\tPASS\t2,1\t25,28\tchr5:70950921-70953320:+\n"));
    let summary = d.collect();
    expect(data(&summary, "parascopy"), "smn1.cn", "1");
    expect(
        data(&summary, "parascopy"),
        "reading",
        r#""one SMN1 copy: the pattern of an SMA carrier""#,
    );
    d.put(SMN, format!("{SMN_HEADER}chr5\t70950920\t70953320\tSMN1\tPASS\t3\t35\tPASS\t1,2\t15,16\tchr5:70075496-70077900:+\n"));
    expect(data(&d.collect(), "parascopy"), "smn1.cn", "null");

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
    expect(find("BRCA2").unwrap(), "weak_call", "true");
    expect(find("LDLR").unwrap(), "weak_call", "true");
    expect(find("MYBPC3").unwrap(), "weak_call", "false");
    expect(find("MYH7").unwrap(), "frequency", r#""unknown""#);

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
    const BLOODTYPE: &str = "bloodtype/S_bloodtype.tsv";
    let d = Sample::new("blood-groups", "S");
    let summary = d.collect();
    assert!(
        matches!(at(&summary, "not_run"), Value::Array(v) if v.iter().any(|x| matches!(x, Value::Str(s) if s == "Blood groups (predicted)")))
    );
    d.put(BLOODTYPE, "system\tphenotype\tgenotype\tconfidence\tnote\nABO\tA1\tA1/O.01\thigh\t\nRhD\tD-\tRHD deleted/RHD deleted\thigh\t\nKell\tK-k+\tk/k\tmedium\tx\n");
    let summary = d.collect();
    expect(section(&summary, "bloodtype"), "step", r#""41""#);
    expect(data(&summary, "bloodtype"), "abo", r#""A1""#);
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
    d.put(
        ROH,
        "RG\tS\tchrY\t2781480\t47781479\t45000000\t900\t50\n\
         RG\tS\tchrEBV\t1\t6000000\t6000000\t90\t50\n\
         RG\tS\tchr1\t125200001\t143200000\t18000000\t300\t50\n\
         RG\tS\tchr9\t42300001\t49300000\t7000000\t200\t50\n\
         RG\tS\tchr1\t120000001\t145000000\t25000000\t500\t50\n\
         RG\tS\tchr2\t10000001\t17000000\t7000000\t900\t50\n",
    );
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
    d.put(ROH, "RG\tS\tchr2\t10000001\t13000000\tinf\t900\t50\n");
    expect(section(&d.collect(), "roh"), "state", r#""unreadable""#);
}

#[test]
fn telomere_counts_require_enough_reads_at_the_observed_depth() {
    const TELOMERE: &str = "telomere/S/S/S_summary.tsv";
    let d = Sample::new("telomere", "S");
    for (count, low) in [(40, "true"), (4000, "false")] {
        d.put(TELOMERE, format!("PID\tsample\ttotal_reads\tread_lengths\trepeat_threshold_set\trepeat_threshold_used\tintratel_reads\tgc_bins_for_correction\ttotal_reads_with_tel_gc\ttel_content\nS\ttumor\t620000000\t150\tn=7\t10\t{count}\t48-52\t150000000\t1.5\n"));
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
fn the_readers_agree_on_clinvar_stars_and_combined_allele_notes() {
    let statuses = [
        ("practice_guideline", 4),
        ("reviewed_by_expert_panel", 3),
        ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
        ("criteria_provided,_single_submitter", 1),
        ("criteria_provided,_conflicting_interpretations", 1),
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
    let Value::Array(hits) = at(data(&summary, "clinvar"), "hits") else {
        panic!("missing hits")
    };
    for hit in hits {
        let Value::Int(pos) = at(hit, "pos") else {
            panic!("non-numeric position")
        };
        assert_eq!(at(hit, "stars"), &Value::Int(statuses[*pos as usize - 1].1));
    }

    let variant = "chr22:30000038:GCCAT:CATGGG";
    let base =
        "GENEINFO=SYNTH1:1;CLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter";
    for (gt, caveat) in [("0/1", true), ("0|1", false), ("1/1", false)] {
        let rows = format!(
            "{VCF}chr22\t30000037\t1\tAGC\tA\t50\tPASS\t{base};PGP_PART_OF={variant}\tGT\t{gt}\n\
            chr22\t30000038\t3\tGCCAT\tCATGGG\t.\tPASS\t{base};PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG\tGT\t{gt}\n"
        );
        d.put(CLINVAR, &rows);
        let mut hit_rows = Vec::new();
        og_clinvar::hits::hits(rows.as_bytes(), &mut hit_rows).unwrap();
        let hit_rows = String::from_utf8(hit_rows).unwrap();
        let summary = d.collect();
        let Value::Array(hits) = at(data(&summary, "clinvar"), "hits") else {
            panic!("missing hits")
        };
        for hit in hits {
            let Value::Str(note) = at(hit, "combined_allele") else {
                panic!("missing combined allele note")
            };
            assert!(!note.is_empty() && hit_rows.contains(note));
            assert_eq!(note.contains("phase"), caveat, "{note}");
        }
    }
}

#[test]
fn hla_concordance_and_kir_keep_the_comparison_and_no_read_reasons() {
    let bash = Sample::new("hla-qc", "S");
    put_hla_qc(&bash, false);
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
    expect(pgx, "low_gq.0.gene", r#""CYP2C19""#);
    expect(pgx, "low_gq.0.min_gq", "12");
    expect(&summary, "cyp2d6.probable.diplotype", r#""*1/*2""#);
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
    expect(&d.collect(), "cyp2d6.probable", "null");
}

#[test]
fn paraphase_counts_other_regions_and_compares_an_smn1_call() {
    let d = Sample::new("paraphase", "S");
    put_paraphase(&d, 1);
    let summary = d.collect();
    let paralogs = data(&summary, "paralogs");
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
    let locus = |name: &str| {
        loci.iter()
            .find(|row| matches!(at(row, "locus"), Value::Str(s) if s == name))
            .unwrap()
    };
    expect(locus("HTT"), "normal_max", r#""26""#);
    expect(locus("HTT"), "zone", r#""27-35""#);
    expect(locus("DMPK"), "zone", r#""35-49""#);
    expect(locus("LOCUS_FL"), "sizing", r#""flanking_only""#);
    expect(locus("LOCUS_SP"), "sizing", r#""sized""#);
    expect(
        locus("LOCUS_LOW"),
        "allele_sizing",
        r#""sized/low_confidence""#,
    );
    expect(locus("RFC1"), "sizing", r#""in_repeat""#);
    expect(expansions, "unsized", r#"["LOCUS_FL","LOCUS_LOW","RFC1"]"#);
    let Value::Array(flagged) = at(expansions, "flagged") else {
        panic!("missing flags")
    };
    let rfc1 = flagged
        .iter()
        .find(|row| matches!(at(row, "locus"), Value::Str(s) if s == "RFC1"))
        .unwrap();
    expect(rfc1, "rfc1_motif.verdict", r#""benign_motif""#);

    d.old(STRANGER);
    // Rewrite ExpansionHunter after the older Stranger output.
    let rows = fs::read(d.path(STR)).unwrap();
    d.put(STR, rows);
    let current = d.collect();
    let expansions = data(&current, "expansions");
    expect(expansions, "stranger", "false");
    expect(expansions, "stranger_outdated", "true");
}

#[test]
fn extra_repeats_keep_the_two_cftr_tracts_and_ugt1a1_agreement() {
    let d = Sample::new("extra-repeats", "S");
    for (five_t, genotype, ugt, agreement) in [
        (true, "TG11/TG12 5T/7T", "*1/*28", "agrees"),
        (false, "TG11-7T/TG11-7T", "*1/*1", "differs"),
    ] {
        put_extra_repeats(&d, five_t);
        let summary = d.collect();
        let expansions = data(&summary, "expansions");
        expect(
            expansions,
            "not_graded",
            r#"{"CFTR":"the CFTR tract","UGT1A1":"the UGT1A1 repeat"}"#,
        );
        let Value::Array(readings) = at(expansions, "readings") else {
            panic!("missing repeat readings")
        };
        let reading = |name: &str| {
            readings
                .iter()
                .find(|row| matches!(at(row, "locus"), Value::Str(s) if s == name))
                .unwrap()
        };
        expect(reading("CFTR"), "five_t", &five_t.to_string());
        expect(reading("CFTR"), "genotype", &format!("\"{genotype}\""));
        expect(reading("UGT1A1"), "genotype", &format!("\"{ugt}\""));
        expect(reading("UGT1A1"), "agreement", &format!("\"{agreement}\""));
        expect(reading("SLC6A4"), "status", r#""not_tested""#);
    }
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

    check(false);
    fs::remove_file(d.path(CONSENSUS_VCF)).unwrap();
    let summary = d.collect();
    expect(data(&summary, "annotsv"), "rows.1.support", "null");
    d.put(CONSENSUS_VCF, b"\x1f\x8b\x08\x00not gzip\n");
    let summary = d.collect();
    expect(data(&summary, "annotsv"), "count", "3");
    assert!(
        matches!(at(data(&summary, "annotsv"), "consensus_note"), Value::Str(s) if s.contains("could not be read"))
    );
    d.put(CHECKED, raw);
    expect(section(&d.collect(), "annotsv"), "state", r#""unreadable""#);
    fs::remove_file(d.path(CHECKED)).unwrap();
    expect(section(&d.collect(), "annotsv"), "state", r#""missing""#);
}
