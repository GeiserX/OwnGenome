//! The report's caveats and the state notes of a section.
//!
//! A caveat (`*_NOTE`) states the limit of the method behind a number. Each
//! card declares the caveats it shows, and the HTML writer prints them under
//! the card's title. The text report prints them in the card's lines. The
//! registry also supplies the same caveats to charts and the live site.

use crate::model::{Malformed, Section, State, py};

pub const CLINVAR_FREQ_NOTE: &str = "Check every P/LP hit's population frequency, a heterozygous hit in a dominant gene as \
    well as a homozygous one: a hit above 1% in gnomAD is unlikely to be a high-penetrance \
    finding, and a hit whose frequency is unknown has not been checked.";
pub const CALLABILITY_NOTE: &str = "A gene listed as not read could not be read reliably: too few of its reads were placed \
    on it with confidence, so a 'no finding' in it says nothing. A gene with a caller of \
    its own is left to that caller's section.";
pub const HLA_NOTE: &str = "HLA typing from short-read WGS is approximate. If either HLA-A or HLA-B allele has a \
    T1K quality of 0 or less, that gene is not passed to PharmCAT (step 36).";
pub const EH_NOTE: &str = "Repeat-expansion calls from short reads can be wrong at some loci, and a long expansion \
    cannot be sized. Disease thresholds: docs/results/structural-and-repeats.md.";
pub const TELOMERE_NOTE: &str = "Relative telomere content (intratelomeric reads per million reads with 48-52% GC), not \
    a telomere length: compare it only with samples of a similar age sequenced on the same \
    platform (docs/10-telomere-analysis.md).";
pub const ANCESTRY_FINE_NOTE: &str = "Each name is a group of reference populations we put together and named, not a \
    nationality. Continental groups are well apart; neighbouring regions are not: splits \
    within Iberia are weak for every method, and Greek and Southern Italian ancestry \
    overlap. Read a regional fraction with its bootstrap range. See \
    docs/26b-ancestry-fine.md.";
pub const TRAITS_NOTE: &str = "Each reading is what the cited studies found on average in the people they studied, not \
    a diagnosis. A row marked ancestry reads differently in other populations; its note \
    says how.";
pub const EXOMISER_NOTE: &str = "A ranking to bring to a clinical geneticist, not a diagnosis: a high score means a \
    rare, likely damaging variant in a gene linked to symptoms like these, not that it \
    causes them, and an empty or low list does not rule a genetic cause out. Coding and \
    splice variants only (docs/48-exomiser.md).";
pub const ANCIENT_NOTE: &str = "Closeness on this PCA is genetic similarity, not descent: a genome near an ancient \
    group shares ancestry with it and does not descend from those people. qpAdm models are \
    not part of step 44; docs/44-ancient.md says how to read the plot.";
pub const COUPLE_NOTE: &str = "A screen, not a diagnosis: the couple report (couple_<a>_<b>.md in the family folder) \
    lists the genes, variant types and repeats it does not see. These counts are of its \
    gene table: a carrier deletion allele from read depth (step 35), and a possible shared \
    risk it adds, is in that report only. See a genetic counsellor before acting on a \
    result, or on its absence (docs/39-couple.md).";
pub const PRS_ANCESTRY_NOTE: &str = "Most PRS scores were developed in European-ancestry populations and predict less well \
    for other ancestries, even with an ancestry-adjusted percentile.";
pub const PHARMCAT_NOTE: &str = "PharmCAT calls star alleles and phenotypes from the variants this genome holds in the \
    genes it covers. A variant it cannot see (a deletion, a duplication, a hybrid gene \
    outside CYP2D6) is not in its answer. Its drug guidance is for a prescriber to read \
    with you: never change a medicine on this report alone.";
pub const CPIC_NOTE: &str = "CPIC recommendations follow the phenotype PharmCAT called. They are guidance for \
    prescribers, not a prescription. A gene that was not called has no recommendation here, \
    and that is not the same as a typical response.";
pub const PGX_CALLERS_NOTE: &str = "CYP2D6 sits beside its look-alike CYP2D7, and short reads often cannot tell them apart, \
    so pypgx, Cyrius and step 36 can disagree or give no call. A disagreement or a missing \
    call means the gene was not read reliably here; a clinical CYP2D6 test settles it.";
pub const CPSR_NOTE: &str = "CPSR grades variants in cancer predisposition genes with ACMG rules. A pathogenic or \
    likely pathogenic row is something to take to a clinical geneticist, not a diagnosis. \
    No row does not rule out an inherited cancer risk: large deletions, genes outside its \
    panel and family history are not part of this screen.";
pub const VNTYPER_NOTE: &str = "The MUC1 repeat is long and GC-rich, and short reads cover it unevenly. VNtyper looks \
    only for known frameshifts in it: a negative result lowers the chance of those and does \
    not rule out ADTKD-MUC1, and a positive one needs a clinical test (SNaPshot for dupC) \
    to confirm it.";
pub const SV_CALLERS_NOTE: &str = "Short reads find large deletions and duplications better than insertions and \
    inversions, and the callers disagree on many events. The totals count every call, most \
    of them common in the population; only the AnnotSV class 4 and 5 rows are graded for \
    disease.";
pub const COVERAGE_NOTE: &str = "Mean depth is an average over the genome: some genes get far fewer reads, and the gene \
    callability card lists those. The sex, identity and contamination checks catch a \
    swapped or mixed sample, not every lab error.";
pub const VARIANTS_NOTE: &str = "Counts of the calls DeepVariant made against the reference. A whole genome usually \
    holds 4 to 5 million, most of them common. PASS calls passed the caller's filter; a \
    count says nothing about health on its own.";
pub const LPA_NOTE: &str = "The KIV-2 repeat count from read depth is a research estimate, not an Lp(a) level: the \
    level is a blood test, and the count explains only part of the difference between \
    people.";
pub const CDNA_SCAN_NOTE: &str = "A locus with reads that jump from one exon to the next is spliced DNA: cDNA or plasmid \
    contamination, or a processed copy of the gene the person carries. Check the reads \
    before trusting a variant call there.";
pub const DEEP_INTRONIC_NOTE: &str = "Known pathogenic deep intronic sites, reported on every run: they lie far from any \
    exon, so most filters miss them. 'ref' with a GQ is a reference call, 'ref (low \
    confidence)' one below the call floor; 'no call' says nothing either way.";
pub const HLA_DERIVED_NOTE: &str = "Information about alleles, not a diagnosis: each line is an association reported in the \
    cited study. Carrying an allele does not mean the condition or reaction will happen, \
    and not carrying it does not rule it out. The genotypes are unphased.";
pub const KIR_NOTE: &str = "KIR genotypes as T1K typed them. No step interprets them: they matter for transplant \
    matching and some immune conditions, not for drug guidance.";
pub const MOSAIC_NOTE: &str = "Loss of chrY in a share of the blood cells of men, and of one chrX in women, becomes \
    more common with age. Population studies link loss of Y with smoking and with higher \
    rates of heart disease and some cancers; loss of X is less studied. That is a \
    population-level association, not a diagnosis for one person. Below about 10% of the \
    cells nothing can be told at 30x (docs/33-sample-qc.md).";
pub const COMPHET_NOTE: &str = "A pair can hit both copies of the gene only in trans. cis: one copy carries both \
    variants. unknown: fewer than 2 reads or read pairs cover both positions, the reads \
    disagree, or there is no BAM; the pairs table gives the reason. A parent's genome or \
    long reads can phase it. Every pair: slivar/<sample>_compound_het_pairs.tsv";
pub const ANNOTSV_NOTE: &str = "A class 4 or 5 SV made by one caller, a deletion on the exon edges of a transcript (the \
    reads of a processed transcript or cDNA aligned across an intron) or one at 1% or more \
    in AnnotSV's benign set is unlikely to be a rare pathogenic SV: check the reads before \
    acting on it (docs/05-annotsv.md). Step 05's rules (assets/sv_rule_regions.tsv and the \
    read-depth check) move a class 4 or 5 call to class 3 when it lies in a coding tandem \
    repeat, an immunoglobulin or T-cell receptor locus or the MHC, names a pseudogene by a \
    disease gene's alias, or its read depth does not show the change; those calls are \
    listed apart with the reason and are not counted here.";
pub const RFC1_NOTE: &str = "RFC1: Stranger's catalog sizes the degenerate AARRG unit and flags a low count, so the \
    motif decides. The motif count reads the reads at the locus; it cannot size a long \
    expansion.";
pub const NOT_GRADED_NOTE: &str = "Not graded: the repeat catalog gives these loci no thresholds (for MUC1, CEL, TYMS and \
    PLIN4 the reference repeat is longer than 300 bases, so no short read spans it), so \
    Stranger gives them no status. Not graded is not normal (docs/09-str-expansions.md).";
pub const PARALOG_NOTE: &str = "Paraphase reads PacBio HiFi reads only. A 'no call' means it saw one haplotype and the \
    depth did not settle one copy or two. A research estimate, not a diagnostic test: \
    confirm a result (for example an HBA deletion or one SMN1 copy) with a clinical \
    laboratory.";
pub const PARASCOPY_NOTE: &str = "Parascopy reads short reads at SMN1 and SMN2 and at the other loci of its catalog. A \
    quality under 20 or a filter other than PASS leaves a copy unresolved. HBA1 and HBA2, \
    CYP21A2 and the gene-dosage loci have no Parascopy model: step 35 counts them from read \
    depth against flanking sequence (CYP21A2 also from six sites where it differs from its \
    pseudogene CYP21A1P), which cannot say which chromosome lost or gained a copy. The \
    copy-number intervals (ClinGen's recurrent CNV regions, carrier deletion alleles, the \
    AZF regions) are read from depth too. A research estimate, not a diagnostic test: \
    confirm a result with a clinical laboratory (MLPA or qPCR; gap-PCR for \
    alpha-thalassaemia).";
pub const PARASCOPY_LOCI_NOTE: &str = "The other loci of step 35's catalog (assets/paralog_loci.tsv): agCN is the copies of \
    the gene and its paralogs together, read only with filter PASS and quality 20 or more; \
    which copy gained or lost one is not said. Not assessed names the reason (most often: \
    no read over the locus).";
pub const ACMG_NOTE: &str = "A list to review with a clinician, not a diagnosis. Each gene's ACMG SF v3.3 reporting \
    rule (assets/acmg_sf.tsv) is applied: one variant in a dominant gene; two in a \
    recessive gene such as MUTYH, BTD or ATP7B, even when their phase is unknown; HFE only \
    as homozygous p.C282Y. A variant the rule does not report is listed apart with the \
    reason. A rare HIGH-impact record is not a ClinVar classification.";
pub const CARRIER_NOTE: &str = "Carrier list: the ClinVar P/LP hits of this sample in the genes of \
    assets/carrier_genes.tsv (autosomal recessive and X-linked conditions), one person, no \
    partner needed. A carrier of one recessive change is usually healthy; the result \
    matters for a child when a partner carries a change in the same gene (step 39). A \
    hypomorph or pseudodeficiency allele is set aside.";
pub const Y_TREE_NOTE: &str = "A second placement, on a public tree release newer than the one Yleaf carries, walked \
    from the reads by og y tree-place. The two can differ in depth and in name: branches \
    get renamed and split between releases, and the position tables differ. A branch is \
    taken when it has more derived than ancestral SNPs, and an ancestral call above it \
    blocks the branches below unless the SNPs below outvote that one site. Every level is \
    in y_haplogroup/<sample>_y_tree_path.tsv; docs/37-y-haplogroup.md says how the walk \
    works.";
pub const YSTR_NOTE: &str = "Repeat counts read from the genome with HipSTR, in FamilyTreeDNA's Y-111 order. A count \
    can differ from FamilyTreeDNA's by a fixed number of repeats at some markers, as labs \
    count some repeat blocks differently: compare with a known result before trusting a \
    difference. Multi-copy markers and repeats longer than the reads are not called. Every \
    marker, with its reads and notes, is in y_str/<sample>_y_str.tsv; \
    y_str/<sample>_y_str.csv is the same list for a spreadsheet.";
pub const BLOODTYPE_NOTE: &str = "Predicted from the DNA at the variants that define each antigen, not typed from blood. \
    RhD comes from the read depth over RHD; weak and partial D are not assessed. The \
    evidence for each site is in bloodtype/<sample>_bloodtype_evidence.tsv.";
pub const ARCHAIC_NOTE: &str = "An estimate from IBDmix against the Altai and Vindija Neanderthals and the Altai \
    Denisovan: a segment counts when either copy carries it, and part of what any person \
    gets is old shared ancestry, not gene flow (Africans get some too). Every segment: \
    archaic/<sample>_archaic_segments.tsv; docs/43-archaic.md says where the numbers end.";
pub const CHIP_NOTE: &str = "Clonal haematopoiesis is common with age and most people with it never get a blood \
    cancer: a call here is not a diagnosis. A call is a change in a CHIP gene's published \
    classes, found by Mutect2 without a matched normal; one at a VAF of 35% or more is left \
    out as possibly inherited. docs/49-chip.md says how to read it.";
pub const LOCITYPER_NOTE: &str = "A second opinion beside the primary caller, never the call itself. Locityper (Prodanov \
    et al. 2025, Nat Genet, PMID 41107551) picks the two haplotypes of the HPRC pangenome \
    panel that best explain the reads at each locus; it names haplotypes of panel people \
    (<sample>.1 or .2), not alleles, so the two callers are shown side by side and not \
    compared. Filter: the authors' HLA and KIR benchmark filter (no warning, weighted \
    distance of 30 or less, fewer than 1000 or 20% unexplained reads). Every column: \
    locityper/<sample>_locityper.tsv.";
pub const VIRAL_NOTE: &str = "From the reads that do not map to the human genome, realigned to a fixed set of human \
    viral genomes. Not a test for an infection: most adults carry herpesviruses that stay \
    in their cells, and reagents carry viral DNA too. HHV-6 integration follows Liu et al. \
    2020 (depth over the viral genome at 0.30 of the human depth or more). Every genome: \
    viral/<sample>_viral.tsv.";
pub const MEI_NOTE: &str = "Alu, L1 and SVA copies the genome carries and the reference does not, found and \
    genotyped by MEGAnE (Kojima et al. 2023) from short read pairs. Most are common and \
    harmless; the catalog rows are insertions known to cause disease. Not callable: too few \
    reads in that window to tell. Not seen is not proof of absence: short reads find an \
    SVA's edges, not its body. Call filtered out: MEGAnE called an insertion there but \
    flagged it (LC is low confidence), so it is not a finding, and not proof of absence \
    either. Every call: mei/<sample>_mei.tsv.";
pub const MCA_NOTE: &str = "Large gains, losses and copy-neutral losses of heterozygosity (CN-LOH) carried by a \
    share of the blood cells, from MoChA on the phased heterozygous sites of chr1 to chr22. \
    In a UK Biobank study they became more common with age and went with a higher risk of a \
    later blood cancer (Loh et al. 2018). That is a population-level association, not a \
    diagnosis; an event is confirmed on a new blood sample. Below the floor nothing is \
    reported, and an event above it can still be missed. Every MoChA call, with the reason \
    it was or was not reported: mosaic/<sample>_mocha.tsv.";
macro_rules! roh_note {
    (group $group:literal) => {
        roh_note!(concat!(
            "Allele frequencies came from the panel's ",
            $group,
            " group."
        ))
    };
    ($af:expr) => {
        concat!(
            "Autosomes only: chrX and chrY are left out, and so is the part of a segment in a \
            centromere or heterochromatin band (GRCh38 cytoband acen, gvar, stalk), where few \
            markers map. ",
            $af
        )
    };
}

pub const ROH_NOTE: &str = roh_note!(
    "No population allele frequencies were used (bcftools roh --AF-dflt 0.4): \
    read the total and the short segments as rough; the long segments are the signal."
);
const ROH_GLOBAL_NOTE: &str = roh_note!("Allele frequencies came from the global panel.");
const ROH_GLOBAL_NO_TABLE_NOTE: &str = roh_note!(
    "Allele frequencies came from the global panel: there was no step 26 ancestry table \
    (no ancestry was run, or PRS failed)."
);
const ROH_GLOBAL_NO_GROUP_NOTE: &str = roh_note!(
    "Allele frequencies came from the global panel: no step 26 ancestry group reached 0.8."
);
const ROH_AFR_NOTE: &str = roh_note!(group "afr");
const ROH_AMR_NOTE: &str = roh_note!(group "amr");
const ROH_EAS_NOTE: &str = roh_note!(group "eas");
const ROH_NFE_NOTE: &str = roh_note!(group "nfe");
const ROH_SAS_NOTE: &str = roh_note!(group "sas");
const ROH_UNKNOWN_AF_NOTE: &str =
    roh_note!("The ROH file names an unrecognized allele-frequency mode.");

/// The ROH caveat identified by the first word of the step's AF header.
/// The flat-AF run has no header, as in summaries from before the panel.
/// A global header says why the global frequencies were used at its end
/// (scripts/11-roh-analysis.sh, modules/local/roh); one without a reason,
/// from an older run, gets the plain global note.
pub fn roh_note_key(af_used: Option<&pycompat::json::Value>) -> &'static str {
    use pycompat::json::Value;
    let Some(Value::Str(af)) = af_used else {
        return "ROH_NOTE";
    };
    match af.split_whitespace().next() {
        Some("global") if af.ends_with("; no step 26 ancestry table") => "ROH_GLOBAL_NO_TABLE_NOTE",
        Some("global") if af.ends_with("; no step 26 group at 0.8 or more") => {
            "ROH_GLOBAL_NO_GROUP_NOTE"
        }
        Some("global") => "ROH_GLOBAL_NOTE",
        Some("afr") => "ROH_AFR_NOTE",
        Some("amr") => "ROH_AMR_NOTE",
        Some("eas") => "ROH_EAS_NOTE",
        Some("nfe") => "ROH_NFE_NOTE",
        Some("sas") => "ROH_SAS_NOTE",
        Some("0.4") => "ROH_NOTE",
        _ => "ROH_UNKNOWN_AF_NOTE",
    }
}
pub const CNVPYTOR_NOTE: &str = "Filtered: e-value below 0.01, under half of the reads with mapping quality 0 (q0) and \
    under half N bases (pN), no overlap with a centromere or heterochromatin band, and no \
    chrY call when the declared sex is female.";
pub const MITO_POLYC_NOTE: &str = "Length heteroplasmy in the D-loop poly-C tracts is found in most people and is usually \
    technical; it is listed here and not counted as heteroplasmy.";
pub const MITO_NOTE: &str = "ClinVar and gnomAD v3.1 (chrM sites) annotate the calls in step 20. MITOMAP \
    (mitomap.org) has the curated disease status of mtDNA variants; it is a manual lookup, \
    as it refuses scripted requests.";
pub const ROH_UPD_NOTE: &str = "Long runs concentrated on one chromosome, with few elsewhere, are the sign of \
    isodisomy: both copies of that chromosome from one parent. They can also be chance in a \
    child of related parents. Which parent gave both copies cannot be told without the \
    parents' samples, and a heterodisomy (two different copies from one parent) leaves no \
    runs at all. The parents' genotypes confirm a UPD. Where an imprinted region is listed, \
    a methylation test of that region tests for its imprinting disorder; on its own it does \
    not prove a UPD (docs/11-roh-analysis.md).";
pub const ROH_F_NOTE: &str =
    "F_ROH is the share of the autosomes (outside the bands) in runs of 5 MB or more.";
pub const ROH_BAND_NOTE: &str = "The band is the row of docs/11-roh-analysis.md's table it is nearest to: an average, \
    and neighbouring rows overlap, so it is consistent with that relationship, not a proof \
    of it.";
pub const KARYOTYPE_NOTE: &str = "Named from the copy number of chrX and chrY in the reads. A mixed sample or a sample \
    swap can give the same numbers. A mosaic in about 30 to 70% of the cells gives numbers \
    between whole ones and no name; one in more of the cells can round to a whole number \
    and be named (a Y lost in 75% of the cells reads 45,X). Confirm with a karyotype or a \
    chromosomal microarray (docs/33-sample-qc.md).";
pub const MITO_DELETION_NOTE: &str = "A single large mtDNA deletion is often found in muscle and not in blood (except in \
    Pearson syndrome), so a blood sample with none does not rule one out \
    (docs/20-mtoolbox.md).";
pub const MITO_COPIES_NOTE: &str = "mtDNA copies per cell: twice the chrM depth over the autosomal depth, from the read \
    counts of the BAM index, for a whole-genome BAM (an exome gives no such number). It \
    depends on the tissue and the library, so it compares samples of one kind only.";
pub const VIEWER_NOTE: &str = "    <p class=\"viewer\">The Viewer column opens the reads at each row in the genome \
    viewer on this computer (start it first: scripts/owngenome-view.sh &lt;sample&gt;, \
    docs/viewer.md).</p>";
pub const EXPLORE_NOTE: &str = "Nothing is uploaded by this pipeline. The docs link lists the files each service reads, \
    which steps 12 and 37 write when they run, and what an upload shares.";
pub const CONSUMER_NOTE: &str = "These files stay on this machine. Uploading one to a service shares your genotypes with \
    it, and with the relatives it matches you to: docs/40-export-consumer.md says what to \
    weigh first.";

/// Report caveats by name, including the ROH allele-frequency modes.
pub const NOTES: [(&str, &str); 65] = [
    ("CLINVAR_FREQ_NOTE", CLINVAR_FREQ_NOTE),
    ("CALLABILITY_NOTE", CALLABILITY_NOTE),
    ("HLA_NOTE", HLA_NOTE),
    ("EH_NOTE", EH_NOTE),
    ("TELOMERE_NOTE", TELOMERE_NOTE),
    ("ANCESTRY_FINE_NOTE", ANCESTRY_FINE_NOTE),
    ("TRAITS_NOTE", TRAITS_NOTE),
    ("EXOMISER_NOTE", EXOMISER_NOTE),
    ("ANCIENT_NOTE", ANCIENT_NOTE),
    ("COUPLE_NOTE", COUPLE_NOTE),
    ("PRS_ANCESTRY_NOTE", PRS_ANCESTRY_NOTE),
    ("PHARMCAT_NOTE", PHARMCAT_NOTE),
    ("CPIC_NOTE", CPIC_NOTE),
    ("PGX_CALLERS_NOTE", PGX_CALLERS_NOTE),
    ("CPSR_NOTE", CPSR_NOTE),
    ("VNTYPER_NOTE", VNTYPER_NOTE),
    ("SV_CALLERS_NOTE", SV_CALLERS_NOTE),
    ("COVERAGE_NOTE", COVERAGE_NOTE),
    ("VARIANTS_NOTE", VARIANTS_NOTE),
    ("LPA_NOTE", LPA_NOTE),
    ("CDNA_SCAN_NOTE", CDNA_SCAN_NOTE),
    ("DEEP_INTRONIC_NOTE", DEEP_INTRONIC_NOTE),
    ("HLA_DERIVED_NOTE", HLA_DERIVED_NOTE),
    ("KIR_NOTE", KIR_NOTE),
    ("MOSAIC_NOTE", MOSAIC_NOTE),
    ("COMPHET_NOTE", COMPHET_NOTE),
    ("ANNOTSV_NOTE", ANNOTSV_NOTE),
    ("RFC1_NOTE", RFC1_NOTE),
    ("NOT_GRADED_NOTE", NOT_GRADED_NOTE),
    ("PARALOG_NOTE", PARALOG_NOTE),
    ("PARASCOPY_NOTE", PARASCOPY_NOTE),
    ("PARASCOPY_LOCI_NOTE", PARASCOPY_LOCI_NOTE),
    ("ACMG_NOTE", ACMG_NOTE),
    ("CARRIER_NOTE", CARRIER_NOTE),
    ("Y_TREE_NOTE", Y_TREE_NOTE),
    ("YSTR_NOTE", YSTR_NOTE),
    ("BLOODTYPE_NOTE", BLOODTYPE_NOTE),
    ("ARCHAIC_NOTE", ARCHAIC_NOTE),
    ("CHIP_NOTE", CHIP_NOTE),
    ("LOCITYPER_NOTE", LOCITYPER_NOTE),
    ("VIRAL_NOTE", VIRAL_NOTE),
    ("MEI_NOTE", MEI_NOTE),
    ("MCA_NOTE", MCA_NOTE),
    ("ROH_NOTE", ROH_NOTE),
    ("ROH_GLOBAL_NOTE", ROH_GLOBAL_NOTE),
    ("ROH_GLOBAL_NO_TABLE_NOTE", ROH_GLOBAL_NO_TABLE_NOTE),
    ("ROH_GLOBAL_NO_GROUP_NOTE", ROH_GLOBAL_NO_GROUP_NOTE),
    ("ROH_AFR_NOTE", ROH_AFR_NOTE),
    ("ROH_AMR_NOTE", ROH_AMR_NOTE),
    ("ROH_EAS_NOTE", ROH_EAS_NOTE),
    ("ROH_NFE_NOTE", ROH_NFE_NOTE),
    ("ROH_SAS_NOTE", ROH_SAS_NOTE),
    ("ROH_UNKNOWN_AF_NOTE", ROH_UNKNOWN_AF_NOTE),
    ("CNVPYTOR_NOTE", CNVPYTOR_NOTE),
    ("MITO_POLYC_NOTE", MITO_POLYC_NOTE),
    ("MITO_NOTE", MITO_NOTE),
    ("ROH_UPD_NOTE", ROH_UPD_NOTE),
    ("ROH_F_NOTE", ROH_F_NOTE),
    ("ROH_BAND_NOTE", ROH_BAND_NOTE),
    ("KARYOTYPE_NOTE", KARYOTYPE_NOTE),
    ("MITO_DELETION_NOTE", MITO_DELETION_NOTE),
    ("MITO_COPIES_NOTE", MITO_COPIES_NOTE),
    ("VIEWER_NOTE", VIEWER_NOTE),
    ("EXPLORE_NOTE", EXPLORE_NOTE),
    ("CONSUMER_NOTE", CONSUMER_NOTE),
];

/// The text of the caveat called `name`, or `None` when there is none.
pub fn text(name: &str) -> Option<&'static str> {
    NOTES.iter().find(|(n, _)| *n == name).map(|(_, t)| *t)
}

/// The caveats of each section, by name, in the order the card shows them:
/// the first is the one the live site shows under a card's title.
pub const SECTION_NOTES: [(&str, &[&str]); 37] = [
    ("qc", &["COVERAGE_NOTE", "KARYOTYPE_NOTE"]),
    ("callability", &["CALLABILITY_NOTE"]),
    ("variants", &["VARIANTS_NOTE"]),
    (
        "clinvar",
        &["CLINVAR_FREQ_NOTE", "CDNA_SCAN_NOTE", "DEEP_INTRONIC_NOTE"],
    ),
    ("clinical-filter", &["ACMG_NOTE"]),
    ("cpsr", &["CPSR_NOTE"]),
    ("slivar", &["COMPHET_NOTE"]),
    ("exomiser", &["EXOMISER_NOTE"]),
    (
        "structural-variants",
        &["SV_CALLERS_NOTE", "ANNOTSV_NOTE", "CNVPYTOR_NOTE"],
    ),
    ("mei", &["MEI_NOTE"]),
    ("expansions", &["EH_NOTE", "NOT_GRADED_NOTE", "RFC1_NOTE"]),
    ("vntyper", &["VNTYPER_NOTE"]),
    (
        "mito",
        &[
            "MITO_NOTE",
            "MITO_POLYC_NOTE",
            "MITO_DELETION_NOTE",
            "MITO_COPIES_NOTE",
        ],
    ),
    ("lpa", &["LPA_NOTE"]),
    ("chip", &["CHIP_NOTE"]),
    ("mosaic", &["MOSAIC_NOTE"]),
    ("mosaic-mca", &["MCA_NOTE"]),
    ("carrier-list", &["CARRIER_NOTE"]),
    ("couple", &["COUPLE_NOTE"]),
    (
        "paralogs",
        &["PARASCOPY_NOTE", "PARALOG_NOTE", "PARASCOPY_LOCI_NOTE"],
    ),
    ("pharmcat", &["PHARMCAT_NOTE"]),
    ("cpic", &["CPIC_NOTE"]),
    ("cyp2d6", &["PGX_CALLERS_NOTE"]),
    ("hla", &["HLA_NOTE", "HLA_DERIVED_NOTE", "KIR_NOTE"]),
    ("locityper", &["LOCITYPER_NOTE"]),
    ("ancestry-fine", &["ANCESTRY_FINE_NOTE"]),
    ("mt-haplogroup", &["EXPLORE_NOTE"]),
    (
        "y-chromosome",
        &["Y_TREE_NOTE", "YSTR_NOTE", "EXPLORE_NOTE"],
    ),
    ("archaic", &["ARCHAIC_NOTE"]),
    ("ancient", &["ANCIENT_NOTE"]),
    (
        "roh",
        &["ROH_NOTE", "ROH_UPD_NOTE", "ROH_F_NOTE", "ROH_BAND_NOTE"],
    ),
    ("traits-panel", &["TRAITS_NOTE"]),
    ("prs", &["PRS_ANCESTRY_NOTE"]),
    ("bloodtype", &["BLOODTYPE_NOTE"]),
    ("telomere", &["TELOMERE_NOTE"]),
    ("viral", &["VIRAL_NOTE"]),
    ("consumer-export", &["CONSUMER_NOTE"]),
];

/// The caveat names of section `sid`, empty for an unknown one.
pub fn section_notes(sid: &str) -> &'static [&'static str] {
    SECTION_NOTES
        .iter()
        .find(|(s, _)| *s == sid)
        .map_or(&[], |(_, n)| n)
}

/// The note of a section whose result is not this run's own: from an
/// earlier run, unreadable or failed; `""` for the others.
pub fn stale_note(sec: &Section<'_>) -> Result<String, Malformed> {
    let lead = match sec.state {
        State::Stale => "STALE: ",
        State::Unreadable => "could not be read: ",
        State::Failed => "FAILED: ",
        _ => return Ok(String::new()),
    };
    Ok(format!("{lead}{}", py(sec.note)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chapters::CHAPTERS;

    #[test]
    fn every_section_has_its_caveats_and_every_name_has_a_text() {
        let sids: Vec<&str> = CHAPTERS
            .iter()
            .flat_map(|c| c.sections.iter().map(|s| s.id))
            .collect();
        let listed: Vec<&str> = SECTION_NOTES.iter().map(|(s, _)| *s).collect();
        assert_eq!(sids, listed);
        for (sid, names) in SECTION_NOTES {
            for n in names {
                assert!(text(n).is_some(), "{sid}: no caveat {n}");
            }
        }
    }

    #[test]
    fn a_global_roh_note_says_why_the_global_frequencies_were_used() {
        use pycompat::json::Value;
        let key = |s: &str| roh_note_key(Some(&Value::Str(s.into())));
        let p = "global (panel roh_af_gnomad_v4.1.tsv.gz)";
        assert_eq!(
            key(&format!("{p}; no step 26 ancestry table")),
            "ROH_GLOBAL_NO_TABLE_NOTE"
        );
        assert_eq!(
            key(&format!("{p}; no step 26 group at 0.8 or more")),
            "ROH_GLOBAL_NO_GROUP_NOTE"
        );
        assert_eq!(key(p), "ROH_GLOBAL_NOTE");
        assert_eq!(
            key("nfe (panel roh_af_gnomad_v4.1.tsv.gz); step 26 gives EUR 0.93"),
            "ROH_NFE_NOTE"
        );
        assert_eq!(roh_note_key(None), "ROH_NOTE");
        for k in ["ROH_GLOBAL_NO_TABLE_NOTE", "ROH_GLOBAL_NO_GROUP_NOTE"] {
            let t = text(k).unwrap();
            assert!(
                t.starts_with(&ROH_GLOBAL_NOTE[..ROH_GLOBAL_NOTE.len() - 1]),
                "{k}: {t}"
            );
        }
    }
}
