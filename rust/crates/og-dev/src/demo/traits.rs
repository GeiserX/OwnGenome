//! The demo sample's traits chapter and the tables beside it: the traits
//! panel (step 42), blood groups (41), the raw-data export (40), viral DNA
//! in the unmapped reads (47) and telomere content (10).

use super::{Ctx, Rng, Scratch, arg, cell, put, read, table};

/// Step 10: TelomereHunter's summary line, with a made-up telomere content.
pub fn telomere(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    put(
        &ctx.path(&format!("telomere/{s}/{s}_summary.tsv")),
        &format!(
            "PID\tsample\ttotal_reads\tread_lengths\trepeat_threshold_set\trepeat_threshold_used\tintratelomeric_reads\tgc_bins_for_correction\ttotal_reads_with_tel_gc\ttel_content\n{s}\ttumor\t825812236\t151\tn=7 per 100 bp\t10\t61873\t48-52\t185204533\t{:.2}\n",
            rng.uniform(260.0, 480.0)
        ),
    )
}

/// Step 40's summary table, in `og consumer export`'s rows: the pinned site
/// table keeps 614,138 of the chip list's 638,463 sites
/// (docs/40-export-consumer.md); the counts are made up.
pub fn consumer_export(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let (sites, from_blocks, from_records) = (614_138, 470_912, 140_206);
    let mut reasons = [
        ("no gVCF record", 1954),
        ("inside a deletion or a longer allele", 611),
        ("filtered call", 309),
        ("no-call in the gVCF", 146),
    ];
    let called = from_blocks + from_records;
    if sites - called != reasons.iter().map(|(_, n)| n).sum::<i64>() {
        return Err("the no-call reasons do not add up to the sites not called".to_owned());
    }
    let s = &ctx.sample;
    let mut rows: Vec<(String, String)> = [
        ("sites", sites.to_string()),
        ("called", called.to_string()),
        ("no_call", (sites - called).to_string()),
        ("from_reference_blocks", from_blocks.to_string()),
        ("from_variant_records", from_records.to_string()),
        ("reference_checked", from_records.to_string()),
        ("reference_mismatch", "0".to_owned()),
        ("template_sites", "638463".to_owned()),
        ("table_dropped", "24325".to_owned()),
        ("sex", ctx.sex().to_owned()),
        ("file_23andme", format!("{s}_23andme_v5.zip")),
        ("file_ancestrydna", format!("{s}_ancestrydna_layout.zip")),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v))
    .collect();
    reasons.sort();
    rows.extend(
        reasons
            .iter()
            .map(|(k, v)| (format!("no_call: {k}"), v.to_string())),
    );
    let mut text = String::from("key\tvalue\n");
    for (k, v) in rows {
        text.push_str(&format!("{k}\t{v}\n"));
    }
    put(
        &ctx.path(&format!("consumer_export/{s}_consumer_export.tsv")),
        &text,
    )
}

/// Step 41's rows (`og bloodtype report`'s columns and spelling), a common
/// European set: system, phenotype, genotype, confidence, note. The
/// HPA-7 to HPA-14 and HPA-16 to HPA-35 rows are all `a+b-` and are added
/// by [`bloodtype_rows`].
const BLOODTYPE: &[(&str, &str, &str, &str, &str)] = &[
    ("ABO", "A1", "A1/O.01", "high", ""),
    (
        "RhD",
        "D+",
        "RHD/RHD (RHD copies by read depth 1.96)",
        "high",
        "weak and partial D (RHD missense alleles) are not assessed",
    ),
    ("RhCE", "C+c+ E-e+", "C/c e/e", "high", ""),
    (
        "Kell",
        "K-k+ Kp(a-b+) Js(a-b+)",
        "k/k Kpb/Kpb Jsb/Jsb",
        "high",
        "",
    ),
    ("Duffy", "Fy(a+b+)", "FY*A/FY*B", "high", ""),
    ("Kidd", "Jk(a+b+)", "Jka/Jkb", "high", ""),
    ("MNS", "M+N+ S-s+", "M/N s/s", "high", ""),
    ("Lutheran", "Lu(a-b+)", "Lub/Lub", "high", ""),
    ("Diego", "Di(a-b+)", "Dib/Dib", "high", ""),
    ("Dombrock", "Do(a+b+)", "Doa/Dob", "medium", ""),
    ("Colton", "Co(a+b-)", "Coa/Coa", "high", ""),
    ("Vel", "Vel+", "Vel/Vel", "high", ""),
    ("Langereis", "Lan+", "Lan/Lan", "high", ""),
    ("Junior", "Jr(a+)", "Jra/Jra", "high", ""),
    ("H", "H+", "H/H", "high", ""),
    (
        "Secretor",
        "secretor",
        "Se/se",
        "high",
        "null allele: FUT2 c.461G>A x1; weak secretor alleles (such as FUT2 c.418A>T, also written c.385A>T) and a deletion of the whole FUT2 gene are not assessed",
    ),
    ("HPA-1", "HPA-1(a+b-)", "1a/1a", "high", ""),
    ("HPA-2", "HPA-2(a+b-)", "2a/2a", "high", ""),
    ("HPA-3", "HPA-3(a+b+)", "3a/3b", "high", ""),
    ("HPA-4", "HPA-4(a+b-)", "4a/4a", "high", ""),
    ("HPA-5", "HPA-5(a+b-)", "5a/5a", "high", ""),
    ("HPA-6", "HPA-6(a+b-)", "6a/6a", "high", ""),
    // HPA-7 to HPA-14 here
    ("HPA-15", "HPA-15(a+b+)", "15a/15b", "high", ""),
    // HPA-16 to HPA-35 here
    (
        "HNA-1",
        "HNA-1(a+b+c-)",
        "FCGR3B*01/FCGR3B*02",
        "low",
        "FCGR3B and FCGR3A are near copies and the FCGR3B copy number varies (none to four or more): short reads can misplace these calls, the copies are not counted, and HNA-1null (no FCGR3B) is not assessed",
    ),
    ("HNA-3", "HNA-3(a+b-)", "3a/3a", "high", ""),
    ("HNA-4", "HNA-4(a+b-)", "4a/4a", "high", ""),
    ("HNA-5", "HNA-5(a+b-)", "5a/5a", "high", ""),
];

/// [`BLOODTYPE`] with the HPA rows it leaves out, in step 41's order:
/// system, phenotype, genotype, confidence, note.
fn bloodtype_rows() -> Vec<[String; 5]> {
    let plain = |n: u32| {
        [
            format!("HPA-{n}"),
            format!("HPA-{n}(a+b-)"),
            format!("{n}a/{n}a"),
            "high".to_owned(),
            String::new(),
        ]
    };
    let mut out = Vec::new();
    for (system, phen, gt, conf, note) in BLOODTYPE {
        out.push([
            (*system).to_owned(),
            (*phen).to_owned(),
            (*gt).to_owned(),
            (*conf).to_owned(),
            (*note).to_owned(),
        ]);
        match *system {
            "HPA-6" => out.extend((7..15).map(plain)),
            "HPA-15" => out.extend((16..36).map(plain)),
            _ => {}
        }
    }
    out
}

/// Step 41: the rows of [`bloodtype_rows`], checked against the systems
/// `og bloodtype systems` names. A platelet or neutrophil antigen the
/// phenotype lacks gets the marker table's meaning at the head of the
/// note, as `og bloodtype report` writes it.
pub fn bloodtype(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let rows = bloodtype_rows();
    let names = ctx.og_output(&["bloodtype", "systems"])?;
    let names: Vec<&str> = names.lines().collect();
    let ours: Vec<&str> = rows.iter().map(|r| r[0].as_str()).collect();
    if ours != names {
        return Err(format!(
            "the demo lists the blood systems {ours:?}, step 41 writes {names:?}"
        ));
    }
    let markers = table(&read(&ctx.asset("assets/bloodtype_markers.tsv"))?);
    let mut means: Vec<(String, String)> = Vec::new();
    for r in &markers {
        let (system, meaning) = (cell(r, "system"), cell(r, "meaning"));
        if !meaning.is_empty() && !means.iter().any(|(s, _)| s == system) {
            means.push((system.to_owned(), meaning.to_owned()));
        }
    }
    let mut text = String::from("system\tphenotype\tgenotype\tconfidence\tnote\n");
    for [system, phen, gt, conf, note] in rows {
        // The antigens the phenotype lacks: the letter before each "-" inside its brackets.
        let inside = phen.split_once('(').map_or("", |(_, r)| r);
        let b = inside.as_bytes();
        let lacks: Vec<char> = (1..b.len())
            .filter(|&i| b[i] == b'-' && b[i - 1].is_ascii_lowercase())
            .map(|i| b[i - 1] as char)
            .collect();
        let mut note = note;
        if let Some((_, meaning)) = means.iter().find(|(s, _)| *s == system)
            && !lacks.is_empty()
        {
            let (prefix, num) = system.split_once('-').unwrap_or((system.as_str(), ""));
            let names: Vec<String> = lacks.iter().map(|a| format!("{prefix}-{num}{a}")).collect();
            let head = format!("lacks {}: {meaning}", names.join(", "));
            note = if note.is_empty() {
                head
            } else {
                format!("{head}; {note}")
            };
        }
        text.push_str(&format!("{system}\t{phen}\t{gt}\t{conf}\t{note}\n"));
    }
    put(
        &ctx.path(&format!("bloodtype/{}_bloodtype.tsv", ctx.sample)),
        &text,
    )
}

/// Step 42: copies of each row's alt allele, a common European set; the
/// TAS2R38 rows are one taster and one non-taster haplotype, APOE is e3/e3.
/// AHR and OR2M7 (rs4410790, rs4481887) have no genotype, as a gap in the
/// gVCF.
const TRAIT_ALT_COPIES: &[(&str, usize)] = &[
    ("rs4988235", 1),
    ("rs182549", 1),
    ("rs145946881", 0),
    ("rs41380347", 0),
    ("rs41525747", 0),
    ("rs869051967", 0),
    ("rs671", 0),
    ("rs1229984", 2),
    ("rs762551", 2),
    ("rs2472297", 1),
    ("rs713598", 1),
    ("rs1726866", 1),
    ("rs10246939", 1),
    ("rs17822931", 0),
    ("rs10427255", 1),
    ("rs1815739", 1),
    ("rs12913832", 1),
    ("rs1800407", 0),
    ("rs16891982", 2),
    ("rs1426654", 0),
    ("rs12203592", 0),
    ("rs1805007", 0),
    ("rs1805008", 0),
    ("rs1805009", 0),
    ("rs3827760", 0),
    ("rs1801133", 1),
    ("rs601338", 1),
    ("rs2282679", 1),
    ("rs887829", 1),
    ("rs1801272", 0),
    ("rs429358", 0),
    ("rs7412", 0),
    ("rs10455872", 0),
    ("rs3798220", 0),
    ("rs78655421", 0),
];

/// Rows of the panel left without a genotype on purpose.
const TRAIT_GAPS: [&str; 2] = ["rs4410790", "rs4481887"];

/// Step 42: [`TRAIT_ALT_COPIES`] at the panel's positions, read by
/// `og traits table`. A panel row the demo gives no genotype and does not
/// leave out on purpose stops the writer, so a new trait gets a value.
pub fn traits(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let panel_path = ctx.asset("assets/traits.tsv");
    let panel = table(&read(&panel_path)?);
    let copies = |rsid: &str| {
        TRAIT_ALT_COPIES
            .iter()
            .find(|(r, _)| *r == rsid)
            .map(|(_, n)| *n)
    };
    let mut unknown: Vec<&str> = panel
        .iter()
        .map(|r| cell(r, "rsid"))
        .filter(|r| copies(r).is_none() && !TRAIT_GAPS.contains(r))
        .collect();
    if !unknown.is_empty() {
        unknown.sort_unstable();
        return Err(format!(
            "assets/traits.tsv has rows the demo gives no genotype: {}",
            unknown.join(", ")
        ));
    }
    let work = Scratch::new("traits")?;
    let genotypes = work.join("g.tsv");
    let mut text = String::new();
    for r in &panel {
        let Some(n) = copies(cell(r, "rsid")) else {
            continue;
        };
        let (rf, alt) = (cell(r, "ref"), cell(r, "alt"));
        let gt: Vec<&str> = std::iter::repeat_n(rf, 2 - n)
            .chain(std::iter::repeat_n(alt, n))
            .collect();
        text.push_str(&format!(
            "{}\t{}\t{rf}\t{alt}\t{}\t60\n",
            cell(r, "chrom"),
            cell(r, "pos"),
            gt.join("/")
        ));
    }
    put(&genotypes, &text)?;
    std::fs::create_dir_all(ctx.path("traits")).map_err(|e| e.to_string())?;
    let rules = ctx.asset("assets/genotype_rules.tsv");
    let out = ctx.path(&format!("traits/{}_traits.tsv", ctx.sample));
    let rules_out = work.join("rules.tsv");
    ctx.og(&[
        "traits",
        "table",
        "--table",
        arg(&panel_path),
        "--rules",
        arg(&rules),
        "--genotypes",
        arg(&genotypes),
        "--sample",
        &ctx.sample,
        "--out",
        arg(&out),
        "--rules-out",
        arg(&rules_out),
    ])
}

/// Step 47: what `og viral report` writes for made-up coverage: no HHV-6,
/// 7 Epstein-Barr reads (under the 10 reads of a call), and the lab's phiX
/// control, at a human depth of 31.
pub fn viral(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let rows: [[&str; 11]; 4] = [
        [
            "HHV-6A",
            "NC_001664.4",
            "159378",
            "0",
            "0.0000",
            "0.000",
            "0",
            "0.000",
            ".",
            "not found",
            "",
        ],
        [
            "HHV-6B",
            "NC_000898.1",
            "162114",
            "0",
            "0.0000",
            "0.000",
            "0",
            "0.000",
            ".",
            "not found",
            "",
        ],
        [
            "Human gammaherpesvirus 4",
            "NC_007605.1",
            "171823",
            "7",
            "0.0041",
            "0.006",
            "0",
            "0.000",
            ".",
            "below threshold",
            "7 counted reads, under the 10 a call needs: too few to tell the virus from contamination at this depth.",
        ],
        [
            "Escherichia phage phiX174",
            "NC_001422.1",
            "5386",
            "1412",
            "1.0000",
            "38.200",
            "37",
            "1.194",
            "2.4645",
            "lab control (phiX)",
            "Illumina adds phiX to the run as a control: these reads were added in the lab.",
        ],
    ];
    let mut text = String::from(
        "virus\taccession\tlength\treads\tbreadth\tmean_depth\tmode_depth\thuman_depth\trelative_depth\tcopies_per_cell\tcall\tnote\n",
    );
    for r in rows {
        text.push_str(&format!(
            "{}\t31.00\t{}\n",
            r[..7].join("\t"),
            r[7..].join("\t")
        ));
    }
    put(&ctx.path(&format!("viral/{}_viral.tsv", ctx.sample)), &text)
}
