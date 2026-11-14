//! Group b of `og summary collect`: the clinical sections. ClinVar, the
//! carrier list, PharmCAT, CPIC, pypgx, Cyrius, CPSR, the clinical filter
//! with its ACMG secondary-findings tier, slivar, and the `cyp2d6` block.
//!
//! Each input is read once: the ClinVar hits VCF feeds the ClinVar
//! section, the ACMG tier, the carrier list and CPSR's call floor; the
//! clinical VCF feeds the record count and the ACMG tier; AnnotSV's table
//! and step 35's copy-number tables feed both two-hit readings. The
//! ClinVar section (with its VEP lookups) and the clinical section (with
//! the clinical VCF) are read on their own threads.

use super::acmg::{self, ClinicalVcf};
use super::carrier::Carrier;
use super::clinical::{Clinical, Slivar};
use super::clinvar::{ClinVar, Floor, HitsVcf};
use super::cpsr::{self, Cpsr};
use super::model::{Section, Spec, Summary};
use super::pgx::{Cpic, Cyp2d6, Cyrius, PharmCat, Pypgx};
use super::scan::FolderIndex;
use super::sections::{Data, SampleQc};
use super::two_hits::{self, DeletionInputs};
use super::{Folder, Input, Read, read_at, with_mtime};

/// The sections `og summary collect --only b` reads, in summary.json's order.
pub const GROUP_B: [Spec; 9] = [
    Spec {
        key: "clinvar",
        title: "ClinVar screen",
        step: "06",
    },
    Spec {
        key: "carrier",
        title: "Carrier list",
        step: "06",
    },
    Spec {
        key: "pharmcat",
        title: "PharmCAT",
        step: "07",
    },
    Spec {
        key: "cpic",
        title: "CPIC drug-gene recommendations",
        step: "27",
    },
    Spec {
        key: "pypgx",
        title: "pypgx",
        step: "32",
    },
    Spec {
        key: "cyrius",
        title: "Cyrius (CYP2D6)",
        step: "21",
    },
    Spec {
        key: "cpsr",
        title: "Cancer predisposition (CPSR)",
        step: "17",
    },
    Spec {
        key: "clinical",
        title: "Clinical variant filter",
        step: "23",
    },
    Spec {
        key: "slivar",
        title: "Variant prioritization (slivar)",
        step: "31",
    },
];

/// Reads the folder into the summary of group b. An error is a run status
/// or a manifest that is there and cannot be read.
pub fn collect(input: &Input) -> Result<Summary, String> {
    let folder = Folder::open(input)?;
    let sections = read(input, &folder);
    let cyp2d6 = cyp2d6(input, &folder.index, &sections);
    let s = pycompat::mtime::escape(&input.sample);
    // The run block is group a's: the declared sex from the command line,
    // the run status, the manifest or step 33's table, the first that has one.
    let qc_declared = folder
        .index
        .first(&[format!("qc/{s}_sample_qc.tsv")])
        .and_then(|rel| {
            let p = crate::join(&input.sample_dir, rel);
            let rows = super::sections::read_tsv(&p).ok()?;
            SampleQc::read(&p, &rows)
                .ok()?
                .declared_sex()
                .map(str::to_owned)
        });
    let declared_sex = super::declared(input, &folder, qc_declared.as_deref());
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0.0, |d| d.as_secs_f64());
    Ok(Summary {
        sample: input.sample.clone(),
        generated_utc: super::time::utc_minute(now).unwrap_or_default(),
        run: folder.run,
        manifest: folder.manifest,
        declared_sex,
        sections,
        cyp2d6,
    })
}

/// Reads the sections of group b, in [`GROUP_B`] order.
pub fn read(input: &Input, folder: &Folder) -> Vec<Section> {
    let d = input.sample_dir.as_str();
    let index = &folder.index;
    let s = pycompat::mtime::escape(&input.sample);
    let raw = input.sample.as_str();
    let file = |rel: &str| -> (String, String) { (crate::join(d, rel), rel.to_owned()) };
    let find = |pats: &[String]| index.first(pats).map(file);
    let one = |p: String| find(&[p]);

    // The shared inputs, each read once.
    let hits = one(format!("clinvar/{s}_clinvar_hits.vcf")).map(|(p, rel)| HitsVcf::read(p, rel));
    let clinical_table = one(format!("clinical/{s}_clinical_summary.tsv"));
    let clinical_vcf_at = one(format!("clinical/{s}_clinical.vcf.gz"));
    let deletions = DeletionInputs {
        annotsv: two_hits::read_annotsv(
            one(format!("annotsv/{s}_sv_class4_5_checked.tsv")),
            one(format!("sv_merged/{s}_sv_consensus.vcf.gz")).map(|(p, _)| p),
        ),
        step35: two_hits::read_step35(d, raw, &input.assets),
    };
    let side = |rel: &str| one(rel.to_owned());

    let (clinvar, clinical, rest) = std::thread::scope(|scope| {
        let clinvar = scope.spawn(|| -> Read {
            let Some(h) = &hits else { return Ok(None) };
            let data = ClinVar::read(h, d, &side, &s)?;
            with_mtime(&h.path, h.rel.clone(), Data::ClinVar(data))
        });
        let clinical = scope.spawn(|| -> Read {
            let at = clinical_table.clone().or_else(|| clinical_vcf_at.clone());
            let Some((p, rel)) = at else { return Ok(None) };
            let vcf = clinical_vcf_at.as_ref().map(|(v, _)| ClinicalVcf::read(v));
            let tier = acmg::Inputs {
                assets: &input.assets,
                hits: hits.as_ref().map(|h| &h.read),
                clinical: vcf.as_ref(),
                deletions: &deletions,
            };
            let data = Clinical::read(
                clinical_table.as_ref().map(|(p, _)| p.as_str()),
                vcf.as_ref(),
                &tier,
            )?;
            with_mtime(&p, rel, Data::Clinical(data))
        });
        let rest = read_rest(input, &s, index, &find, &hits, &deletions);
        let join = |h: std::thread::ScopedJoinHandle<'_, Read>| {
            h.join()
                .unwrap_or_else(|_| Err("internal error: the section reader stopped".into()))
        };
        (join(clinvar), join(clinical), rest)
    });
    let [carrier, pharmcat, cpic, pypgx, cyrius, cpsr, slivar] = rest;
    let reads = [
        clinvar, carrier, pharmcat, cpic, pypgx, cyrius, cpsr, clinical, slivar,
    ];
    GROUP_B
        .iter()
        .zip(reads)
        .map(|(spec, r)| Section::new(*spec, r, folder.run.as_ref()))
        .collect()
}

/// The `cyp2d6` block: CPIC, pypgx and Cyrius side by side, with step 36's
/// row. `None` when `sections` lack one of the three.
pub fn cyp2d6(input: &Input, index: &FolderIndex, sections: &[Section]) -> Option<Cyp2d6> {
    let s = pycompat::mtime::escape(&input.sample);
    let d = input.sample_dir.as_str();
    let sec = |k: &str| sections.iter().find(|x| x.spec.key == k);
    let consensus = index
        .first(&[format!("pgx_consensus/{s}_pgx_consensus.tsv")])
        .map(|rel| (crate::join(d, rel), rel.to_owned()));
    match (sec("cpic"), sec("pypgx"), sec("cyrius")) {
        (Some(c), Some(p), Some(y)) => Some(Cyp2d6::new(c, p, y, consensus, d)),
        _ => None,
    }
}

/// Finds a section's file: the full path and the path in the folder.
type Find<'a> = dyn Fn(&[String]) -> Option<(String, String)> + 'a;

/// The seven sections read on the calling thread, in order.
fn read_rest(
    input: &Input,
    s: &str,
    index: &FolderIndex,
    find: &Find,
    hits: &Option<HitsVcf>,
    deletions: &DeletionInputs,
) -> [Read; 7] {
    let one = |p: String| find(&[p]);
    let carrier = match hits {
        None => Ok(None),
        Some(h) => Carrier::read(&input.assets, &h.read, deletions)
            .and_then(|c| with_mtime(&h.path, h.rel.clone(), Data::Carrier(c))),
    };
    let pharmcat = pharmcat_at(index, input, s)
        .and_then(|at| read_at(at, |p| PharmCat::read(p).map(Data::PharmCat)));
    let cpic = read_at(one(format!("cpic/{s}_phenotypes.tsv")), |p| {
        Cpic::read(p, &input.sample).map(Data::Cpic)
    });
    let pypgx = read_at(one(format!("pypgx/{s}_pypgx_summary.tsv")), |p| {
        Pypgx::read(p, &input.sample).map(Data::Pypgx)
    });
    let cyrius = read_at(
        find(&[
            format!("cyrius/{s}_cyp2d6.tsv"),
            "cyrius/*_cyp2d6.tsv".to_owned(),
        ]),
        |p| Cyrius::read(p).map(Data::Cyrius),
    );
    let html = find(&[
        format!("cpsr/{s}.cpsr.grch38.html"),
        "cpsr/*.cpsr.*.html".to_owned(),
    ]);
    let tsv = find(&[
        format!("cpsr/{s}.cpsr.grch38.classification.tsv.gz"),
        "cpsr/*.cpsr.*.classification.tsv.gz".to_owned(),
    ]);
    let cpsr = read_at(tsv.clone().or_else(|| html.clone()), |_| {
        let floor = match hits {
            None => Ok(Floor::default()),
            Some(h) => h.read.as_ref().map(|(f, _)| *f).map_err(Clone::clone),
        };
        Cpsr::read(&cpsr::Inputs {
            html: html.as_ref().map(|(p, _)| p.clone()),
            tsv: tsv.as_ref().map(|(p, _)| p.clone()),
            floor,
            cdna: one(format!("cpsr/{s}_cpsr_cdna_check.tsv")).map(|(p, _)| p),
        })
        .map(Data::Cpsr)
    });
    let slivar_table = one(format!("slivar/{s}_slivar_summary.tsv"));
    let slivar_vcf = if slivar_table.is_none() {
        one(format!("slivar/{s}_prioritized.vcf.gz"))
    } else {
        None
    };
    let slivar = read_at(slivar_table.clone().or(slivar_vcf.clone()), |_| {
        Slivar::read(
            slivar_table.as_ref().map(|(p, _)| p.as_str()),
            slivar_vcf.as_ref().map(|(p, _)| p.as_str()),
            &input.sample,
        )
        .map(Data::Slivar)
    });
    [carrier, pharmcat, cpic, pypgx, cyrius, cpsr, slivar]
}

/// PharmCAT's report: the newest of the files its four names match.
fn pharmcat_at(
    index: &FolderIndex,
    input: &Input,
    s: &str,
) -> Result<Option<(String, String)>, String> {
    let d = input.sample_dir.as_str();
    let found: Vec<&str> = [
        format!("pharmcat/{s}.report.json"),
        "pharmcat/*.report.json".to_owned(),
        format!("vcf/{s}.report.json"),
        "vcf/*.report.json".to_owned(),
    ]
    .iter()
    .flat_map(|p| index.all(p))
    .collect();
    let paths: Vec<(String, &str)> = found
        .iter()
        .map(|rel| (crate::join(d, rel), *rel))
        .collect();
    let newest = pycompat::mtime::newest(paths.iter().map(|(p, _)| p.as_str()))
        .map_err(|e| e.to_string())?;
    Ok(newest.map(|p| {
        let rel = paths.iter().find(|(q, _)| q == p).map_or("", |(_, r)| r);
        (p.to_owned(), rel.to_owned())
    }))
}
