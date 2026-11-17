//! The opt-in step sections in both reports: CHIP, viral DNA, MEI, MoChA,
//! VNtyper, Locityper, Exomiser, telomere content and runs of homozygosity,
//! plus the SV caller and cDNA scan blocks other sections compose. Moved
//! from the report checks of
//! tests/test_{chip,viral,mei,mocha,vntyper,locityper,exomiser}_report.py;
//! the byte-for-byte match with render_report.py is the `m-optin-*`
//! report-render parity cases (tests/parity/cases.tsv).
//!
//! The section data is what `og summary collect` writes for the tool tables
//! of each step (og-summary's `tests/optin_summary.rs`).

use og_render::model::Summary;
use og_render::notes::{CHIP_NOTE, MCA_NOTE, MEI_NOTE, VIRAL_NOTE};
use og_render::sections::{cdna_scan, chip, exomiser, mei, mocha, sv_callers, viral};
use og_render::{html, text};

/// A summary with the base sections and `sections`, each
/// `(key, title, step, state, data)` with `data` a JSON object.
fn summary(sections: &[(&str, &str, &str, &str, &str)], extra: &str) -> String {
    let mut s = String::new();
    for (key, title, step, state, data) in sections {
        s += &format!(
            r#", "{key}": {{"title": "{title}", "state": "{state}", "step": "{step}", "source": "{key}/S1_{key}.tsv",
             "file_date": "2026-11-06 11:00:00", "note": null, "data": {data}}}"#
        );
    }
    format!(
        r#"{{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-06 12:00 UTC",
 "run": {{"started_utc": "2026-11-06 10:00 UTC", "declared_sex": "male"}},
 "databases": {{"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null}},
 "manifest": [],
 "sections": {{
  "coverage": {{"title": "Coverage (mosdepth)", "state": "ok", "step": "16b", "source": "qc/S1.mosdepth.summary.txt",
               "file_date": "2026-11-06 11:00:00", "note": null, "data": {{"mean_depth": 31.25}}}}{s}
 }},
 "cyp2d6": null, "not_run": [], "not_assessed": ["MUC1 VNTR frameshifts (ADTKD-MUC1): step 52 (VNtyper, opt-in) gave no reading; no standard caller sees them"]{extra}
}}"#
    )
}

fn parse(json: &str) -> Summary {
    Summary::parse(json).expect("a summary")
}

fn both(json: &str) -> (String, String) {
    let s = parse(json);
    (
        text::report(&s).expect("the text report"),
        html::report(&s).expect("the HTML report"),
    )
}

/// The frame of section `sid` in the HTML report, empty when it has none.
fn frame_of<'a>(page: &'a str, sid: &str) -> &'a str {
    let open = format!("id=\"{sid}\">");
    page.split_once(&open).map_or("", |(_, rest)| {
        rest.split("</section>").next().unwrap_or("")
    })
}

// --- step 49, CHIP ------------------------------------------------------------------

const CHIP_CALL: &str = r#"{"calls": [{"gene": "DNMT3A", "protein_change": "R882H", "vaf": "0.130", "alt_reads": "6",
  "depth": "46", "note": "hotspot R882", "status": "chip"}], "other": 0, "median_depth": 30.0,
  "median_limit": 0.137, "not_covered": ["JAK2"]}"#;

fn chip_summary(data: &str) -> String {
    summary(
        &[("chip", "Clonal haematopoiesis (CHIP)", "49", "ok", data)],
        "",
    )
}

/// The planted check of this package: the limit sentence, word for word.
#[test]
fn the_chip_limit_sentence() {
    let s = parse(&chip_summary(CHIP_CALL));
    assert_eq!(
        chip::limit_line(&s.section("chip")).unwrap(),
        "Detection limit: a clone of VAF about 14% or more is found 8 times in 10 at this depth \
         (median 30x over the genes' coding exons). CHIP starts at a VAF of 2%: smaller clones are \
         likely missed, so no call does not rule one out."
    );
    let none = parse(&chip_summary(
        r#"{"calls": [], "median_depth": null, "median_limit": null, "not_covered": []}"#,
    ));
    assert_eq!(
        chip::limit_line(&none.section("chip")).unwrap(),
        "No CHIP gene has the depth to find a clone: no result either way."
    );
}

#[test]
fn the_chip_lines_show_the_call_the_limit_and_the_uncovered_gene() {
    let s = parse(&chip_summary(CHIP_CALL));
    let lines = chip::lines(&s.section("chip"), "S1").unwrap().join("\n");
    for want in [
        "DNMT3A R882H: VAF 0.130 (6 of 46 reads)",
        "about 14% or more",
        "JAK2",
    ] {
        assert!(lines.contains(want), "{want} not in {lines}");
    }
    let s = parse(&chip_summary(r#"{"calls": [], "median_limit": null}"#));
    assert!(
        chip::lines(&s.section("chip"), "S1")
            .unwrap()
            .iter()
            .any(|l| l == "No CHIP call.")
    );
    let (txt, page) = both(&chip_summary(CHIP_CALL));
    assert!(txt.contains("## Clonal Haematopoiesis (CHIP, step 49)") && txt.contains(CHIP_NOTE));
    assert!(frame_of(&page, "chip").contains("<h2>Clonal Haematopoiesis (CHIP)</h2>"));
}

// --- step 47, viral -----------------------------------------------------------------

fn viral_summary(data: &str) -> String {
    summary(&[("viral", "Viral DNA (step 47)", "47", "ok", data)], "")
}

#[test]
fn the_viral_lines() {
    let data = r#"{"hhv6": "inherited integrated HHV-6B likely", "human_depth": "30.00", "genomes": 5, "found": [
      {"virus": "HHV-6B", "call": "inherited integrated HHV-6 likely", "reads": "400", "breadth": "1.0000",
       "relative_depth": "0.500", "copies_per_cell": "1.0000"}]}"#;
    let s = parse(&viral_summary(data));
    let lines = viral::lines(&s.section("viral")).unwrap();
    assert_eq!(lines[0], "HHV-6: inherited integrated HHV-6B likely");
    assert_eq!(lines.last().unwrap(), VIRAL_NOTE);
    let (_, page) = both(&viral_summary(data));
    assert!(frame_of(&page, "viral").contains("HHV-6B"));

    // No human depth: the line leaves out copies per cell.
    let nohuman = r#"{"hhv6": "not assessed (no human depth)", "human_depth": ".", "genomes": 5, "found": [
      {"virus": "Human gammaherpesvirus 4", "call": "viral DNA found", "reads": "30", "breadth": "1.0000",
       "relative_depth": ".", "copies_per_cell": "."}]}"#;
    let s = parse(&viral_summary(nohuman));
    assert_eq!(
        viral::lines(&s.section("viral")).unwrap()[1],
        "Human gammaherpesvirus 4: viral DNA found (30 reads, 100% of the genome)"
    );
    // A genome under ten reads is listed, never "No read".
    let below = r#"{"hhv6": "no integrated HHV-6", "human_depth": "30.00", "genomes": 5, "found": [
      {"virus": "Hepatitis C virus genotype 1", "call": "below threshold", "reads": "5", "breadth": "0.0544",
       "relative_depth": "0.000", "copies_per_cell": "0.0000"}]}"#;
    let s = parse(&viral_summary(below));
    assert!(
        !viral::lines(&s.section("viral"))
            .unwrap()
            .iter()
            .any(|l| l.starts_with("No read"))
    );
}

/// Follow-up 454: a BAM with no unplaced unmapped read had nothing to look
/// in, so the report says the genomes were not assessed, never that none
/// had a read, in both reports.
#[test]
fn no_unplaced_reads_is_not_assessed_not_no_read() {
    let na = r#"{"hhv6": "not assessed (the BAM has no unplaced unmapped reads)", "human_depth": "30.00",
      "genomes": 24, "found": []}"#;
    let (txt, page) = both(&viral_summary(na));
    assert!(
        txt.contains(
            "  Not assessed on any of the 24 viral genomes: the BAM has no unplaced unmapped reads to look in."
        ),
        "{txt}"
    );
    assert!(!txt.contains("No read on any"), "{txt}");
    let card = frame_of(&page, "viral");
    assert!(
        card.contains(
            "Not assessed on any viral genome: the BAM has no unplaced unmapped reads to look in."
        ),
        "{card}"
    );
    assert!(!card.contains("No read on any"), "{card}");
    // Control: with reads to look in and none found, it is "No read".
    let none =
        r#"{"hhv6": "no integrated HHV-6", "human_depth": "30.00", "genomes": 24, "found": []}"#;
    let (txt, page) = both(&viral_summary(none));
    assert!(txt.contains("No read on any of the 24 viral genomes."));
    assert!(frame_of(&page, "viral").contains("No read on any viral genome."));
}

// --- step 50, MEI -------------------------------------------------------------------

const MEI_DATA: &str = r#"{"calls": 2, "by_element": {"Alu": 1, "SVA": 1}, "sites": [
  {"gene": "FKTN", "element": "SVA", "condition": "Fukuyama congenital muscular dystrophy",
   "inheritance": "autosomal recessive", "status": "seen", "mean_depth": "30.0",
   "call": "chr9:105639660 SVA_F 1 copy PASS"},
  {"gene": "TAF1", "element": "SVA", "condition": "X-linked dystonia-parkinsonism", "inheritance": "X-linked",
   "status": "not seen", "mean_depth": "30.0", "call": "."}]}"#;

#[test]
fn the_mei_lines_and_card() {
    let json = summary(
        &[(
            "mei",
            "Mobile element insertions (MEGAnE)",
            "50",
            "ok",
            MEI_DATA,
        )],
        "",
    );
    let s = parse(&json);
    let lines = mei::lines(&s.section("mei")).unwrap();
    assert_eq!(lines[0], "Insertions called: 2 (Alu 1, SVA 1)");
    assert_eq!(
        lines[1],
        "FKTN SVA (Fukuyama congenital muscular dystrophy; autosomal recessive): seen (depth 30.0, chr9:105639660 SVA_F 1 copy PASS)"
    );
    assert_eq!(lines.last().unwrap(), MEI_NOTE);
    let (_, page) = both(&json);
    let card = frame_of(&page, "mei");
    assert!(page.contains("Mobile Element Insertions") && card.contains("SVA_F"));
    assert!(card.contains("<td>FKTN</td>"), "{card}");
}

// --- step 51, MoChA -----------------------------------------------------------------

fn mocha_data(status: &str, message: &str, calls: &str) -> String {
    format!(
        r#"{{"status": "{status}", "message": "{message}", "depth": "30.0", "het_sites": "18000",
  "het_per_10mb": "3000", "chromosomes": "chr20", "floor_read": "0.020", "floor_min": "0.05", "floor": "0.050",
  "calls_mocha": "2", "calls_reported": "1", "phasing": "input", "not_reported": 1, "calls": [{calls}]}}"#
    )
}

const CNLOH: &str = r#"{"chrom": "chr20", "start": "60000", "end": "20000000", "length_mb": "19.94",
  "arms": "20p, to the p telomere", "type": "CN-LOH", "cell_fraction": "0.150", "n_hets": "6600",
  "lod_baf_phase": "900.0", "reported": "yes", "reason": "."}"#;

fn mosaic(data: &str) -> Vec<String> {
    let json = summary(
        &[(
            "mosaic",
            "Mosaic chromosomal alterations (MoChA)",
            "51",
            "ok",
            data,
        )],
        "",
    );
    mocha::lines(&parse(&json).section("mosaic")).unwrap()
}

#[test]
fn the_mocha_lines() {
    let lines = mosaic(&mocha_data("assessed", ".", CNLOH));
    assert!(
        lines[0].starts_with("Reported at a cell fraction of 0.050 or more at 30.0x"),
        "{}",
        lines[0]
    );
    assert_eq!(
        lines[1],
        "CN-LOH on 20p, to the p telomere, chr20:60000-20000000 (19.94 Mb): about 0.150 of the cells"
    );
    assert_eq!(lines.last().unwrap(), MCA_NOTE);
    let na = mosaic(&mocha_data(
        "not assessed",
        "phased BAF autocorrelation 0.0500 over 0.03",
        "",
    ));
    assert!(na[0].starts_with("Not assessed: phased BAF"), "{}", na[0]);
    let none = mosaic(&mocha_data("assessed", ".", ""));
    assert!(
        none[1].starts_with("No mosaic alteration reported"),
        "{}",
        none[1]
    );
}

// --- step 52, VNtyper ---------------------------------------------------------------

const MUC1: &str = "MUC1 VNTR frameshifts (ADTKD-MUC1): step 52 (VNtyper, opt-in) gave no reading; no standard caller sees them";

fn vntyper(status: &str) -> String {
    let data = format!(
        r#"{{"status": "{status}", "confidence": "High_Precision*", "variant": "59dupC", "flag": "Not flagged",
  "alt_depth": "416", "region_depth": "7110", "depth_score": "0.0585", "vntr_mean_depth": "1075.32",
  "vntr_uncovered_pct": "0.00", "coverage_qc": "PASS",
  "sensitivity": "At 1075x the 100x-or-more band: 160 of 161 runs (99%).", "message": "m", "advice": "a"}}"#
    );
    summary(&[("vntyper", "MUC1 VNTR (VNtyper)", "52", "ok", &data)], "")
}

#[test]
fn the_vntyper_section_in_both_reports() {
    let (txt, page) = both(&vntyper("positive"));
    for want in [
        "## MUC1 VNTR Frameshifts (ADTKD-MUC1, VNtyper)",
        "Result: positive",
        "Variant: 59dupC (Kestrel confidence High_Precision*",
        "160 of 161 runs",
    ] {
        assert!(txt.contains(want), "{want} not in {txt}");
    }
    let card = frame_of(&page, "vntyper");
    assert!(
        card.contains("<h2>MUC1 VNTR Frameshifts (VNtyper)</h2>")
            && card.contains("badge-red\">positive<")
            && card.contains("160 of 161 runs"),
        "{card}"
    );
    let (_, page) = both(&vntyper("negative"));
    assert!(frame_of(&page, "vntyper").contains("badge-green\">negative<"));
    // Without step 52: no card, no block, and both reports list it.
    let (txt, page) = both(&summary(&[], ""));
    assert!(!page.contains("MUC1 VNTR Frameshifts") && !txt.contains("## MUC1 VNTR"));
    assert!(txt.contains(MUC1) && page.contains(MUC1));
}

// --- step 53, Locityper -------------------------------------------------------------

const LOCI: &str = r#"{"called": 2, "passed": 1, "loci": [
  {"locus": "HLA-DQB1", "region": "chr6:32652604-32669077", "status": "called", "genotype": "HG00733.1,NA19240.2",
   "haplotype_1": "HG00733.1", "haplotype_2": "NA19240.2", "quality": "27.8", "total_reads": "900",
   "unexplained_reads": "12", "weighted_distance": "50.000", "warnings": ".", "filter": "fail",
   "reason": "weighted distance 50.00 over 30", "primary": "T1K (step 08): not run"},
  {"locus": "CYP2D6", "region": "chr22:42122285-42152960", "status": "called", "genotype": "a,b",
   "haplotype_1": "HG00733.1", "haplotype_2": "NA19240.2", "quality": "27.8", "total_reads": "900",
   "unexplained_reads": "12", "weighted_distance": "0.200", "warnings": ".", "filter": "pass", "reason": ".",
   "primary": "pypgx (step 32): not run"},
  {"locus": "FCGR", "region": "chr1:161506762-161681047", "status": "not called", "genotype": ".",
   "haplotype_1": ".", "haplotype_2": ".", "quality": ".", "total_reads": ".", "unexplained_reads": ".",
   "weighted_distance": ".", "warnings": ".", "filter": ".", "reason": "Locityper wrote no result for this locus",
   "primary": "no other step calls FCGR"}]}"#;

#[test]
fn the_locityper_section_in_both_reports() {
    let (txt, page) = both(&summary(
        &[("locityper", "Locityper (step 53)", "53", "ok", LOCI)],
        "",
    ));
    for want in [
        "## Locityper Second Opinion at Hard Loci (step 53)",
        "HG00733.1 + NA19240.2",
        "filter fail: weighted distance 50.00",
        "primary:   no other step calls FCGR",
        "not called: Locityper wrote no result",
    ] {
        assert!(txt.contains(want), "{want} not in {txt}");
    }
    let card = frame_of(&page, "locityper");
    for want in [
        "Locityper Second Opinion at Hard Loci</h2>",
        "badge-green\">pass<",
        "badge-yellow\">fail<",
        "<th>Primary caller</th>",
    ] {
        assert!(card.contains(want), "{want} not in {card}");
    }
    let (txt, page) = both(&summary(&[], ""));
    assert!(
        !page.contains("Locityper Second Opinion") && !txt.contains("Locityper Second Opinion")
    );
}

// --- step 48, Exomiser --------------------------------------------------------------

fn variant(v: &str, gt: &str, acmg: &str) -> String {
    format!(
        r#"{{"variant": "{v}", "genotype": "{gt}", "consequence": "missense_variant", "hgvs": "p.(Glu563Ala)",
  "acmg": "{acmg}", "disease": "."}}"#
    )
}

fn gene(rank: &str, name: &str, variants: &[String]) -> String {
    format!(
        r#"{{"rank": "{rank}", "gene": "{name}", "moi": "AD", "exomiser_score": "0.998", "phenotype_score": "0.800",
  "variants": [{}]}}"#,
        variants.join(", ")
    )
}

#[test]
fn the_exomiser_section_in_both_reports() {
    let data = format!(
        r#"{{"terms": [{{"id": "HP:0001156", "label": "Brachydactyly"}}], "genes": [{}, {}, {}]}}"#,
        gene(
            "1",
            "GENEA",
            &[variant("chr10:123256215 T>G", "0/1", "LIKELY_PATHOGENIC")]
        ),
        gene(
            "2",
            "GENEB",
            &[
                variant("chr6:1000 G>A", "0/1", "UNCERTAIN_SIGNIFICANCE"),
                variant("chr6:2000 C>T", "1|0", "PATHOGENIC")
            ]
        ),
        gene(
            "3",
            "GENEC",
            &[variant("chrX:5000 A>G", "1", "UNCERTAIN_SIGNIFICANCE")]
        )
    );
    let json = summary(
        &[(
            "exomiser",
            "Symptom-driven ranking (Exomiser)",
            "48",
            "ok",
            &data,
        )],
        "",
    );
    let (txt, page) = both(&json);
    for want in [
        "1. GENEA (AD) score 0.998",
        "chr6:2000 C>T",
        "not a diagnosis",
    ] {
        assert!(txt.contains(want), "{want} not in {txt}");
    }
    let card = frame_of(&page, "exomiser");
    assert_eq!(card.matches("<tr>").count(), 1 + 4, "{card}");
    assert!(card.contains("Brachydactyly (HP:0001156)") && card.contains("not a diagnosis"));
    // The Viewer cell opens every position of the first chromosome a
    // variant names.
    assert_eq!(
        exomiser::locus_span("chr6:2000 C>T + chr6:900 A>C"),
        Some(("chr6".to_owned(), 900, 2000))
    );
}

// --- telomere content and runs of homozygosity --------------------------------------

#[test]
fn telomere_and_roh_cards_show_without_their_steps() {
    let (txt, page) = both(&summary(&[], ""));
    assert!(frame_of(&page, "telomere").contains("Telomere content (relative)"));
    assert!(frame_of(&page, "roh").contains("Runs of Homozygosity"));
    assert!(!txt.contains("## Telomere content") && !txt.contains("## Runs of Homozygosity"));
    let tel = r#"{"tel_content": "1.5", "intratelomeric_reads": 1234567, "depth": 30.0, "read_floor": 2000,
      "low_reads": true}"#;
    let (txt, page) = both(&summary(
        &[("telomere", "Telomere content", "10", "ok", tel)],
        "",
    ));
    assert!(
        txt.contains("  Intratelomeric reads: 1,234,567 (floor 2,000 at 30.0x)"),
        "{txt}"
    );
    assert!(
        frame_of(&page, "telomere").contains("badge-yellow\">too few telomeric reads"),
        "{page}"
    );
}

// --- the blocks other sections compose ----------------------------------------------

#[test]
fn the_sv_caller_blocks() {
    let cn = r#"{"total": 10, "deletions": 6, "duplications": 4, "significant": 3, "significant_filtered": 1,
      "dropped": {"q0": 1, "gap": 0, "band": 1, "chry_female": 0}, "declared_sex": null}"#;
    let json = summary(
        &[
            (
                "manta",
                "Manta",
                "04",
                "ok",
                r#"{"total": 100, "pass": 80}"#,
            ),
            ("cnvpytor", "CNVpytor", "18", "ok", cn),
        ],
        "",
    );
    let s = parse(&json);
    let w = sv_callers::text_blocks(&s).unwrap();
    assert!(
        w.contains(&"  Total SVs: 100  PASS: 80".to_owned()),
        "{w:?}"
    );
    assert!(
        w.contains(
            &"  Significant after filters: 1 (left out: 1 q0 0.5 or more, 1 over a centromere or heterochromatin band); no declared sex, so chrY calls are kept"
                .to_owned()
        ),
        "{w:?}"
    );
    let c = sv_callers::card(&s).unwrap();
    assert_eq!(c.title, "Structural Variants");
    assert_eq!(
        c.caveats.len(),
        2,
        "the CNVpytor note, then the callers' note"
    );
}

#[test]
fn the_cdna_scan_blocks() {
    let data = r#"{"cdna_scan": {"note": "", "loci": [
      {"chrom": "chr7", "start": 100, "end": 900, "junctions": 3, "reads": 12}]}}"#;
    let json = summary(&[("clinvar", "ClinVar", "06", "ok", data)], "");
    let s = parse(&json);
    let d = s.section("clinvar");
    assert_eq!(
        cdna_scan::lines(&d).unwrap(),
        [
            "cDNA scan of every exon: 1 loci with reads across two or more exon junctions",
            "  chr7:100-900: 3 junctions, 12 reads"
        ]
    );
    assert!(cdna_scan::stat(&d).unwrap().is_some());
    let (lines, caveat) = cdna_scan::html_lines(&d).unwrap();
    assert_eq!(lines.len(), 2);
    assert!(caveat.is_some());
    let not_run = r#"{"cdna_scan": {"note": "not run (no BAM)", "loci": []}}"#;
    let s = parse(&summary(&[("clinvar", "ClinVar", "06", "ok", not_run)], ""));
    let d = s.section("clinvar");
    assert_eq!(
        cdna_scan::lines(&d).unwrap(),
        ["cDNA scan of every exon: not run (no BAM)"]
    );
    assert!(cdna_scan::stat(&d).unwrap().is_none());
}
