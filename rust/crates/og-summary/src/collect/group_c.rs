//! Group c of `og summary collect`: copy number, HLA, polygenic scores,
//! structural variants, repeats, telomeres, runs of homozygosity,
//! haplogroups, the step tables of the opt-in steps, chrM, traits with the
//! APOL1 G2 deletion, and the couple screen.
//!
//! Each section is read on its own, and an error stays in its section. The
//! three parts of the Python collector's post pass that change these
//! sections' data run here too, after every section is read: the repeat
//! readings with the traits and CPIC states, Locityper's primary calls, and
//! the chrY filter of CNVpytor's count for a declared female. What they
//! need from group b (the CPIC state, pypgx and Cyrius) comes in
//! [`Others`]; `--only c` reads none of it, so it is empty there. A part
//! that fails makes only its own section unreadable.

mod chrm;
mod cn;
mod hla;
mod repeats;
pub mod roh;
mod sv;
mod tables;
mod traits;
pub mod util;

use super::model::{Spec, State};
use super::scan::FolderIndex;
use super::sections::Data;
use pycompat::json::Value;

/// The sections of group c, in summary.json's order.
pub const GROUP_C: [Spec; 32] = [
    spec("couple", "Couple carrier screen", "39"),
    spec(
        "parascopy",
        "SMN1/SMN2 copy number (Parascopy, short reads)",
        "35",
    ),
    spec("paralogs", "Paralog genes (Paraphase, PacBio HiFi)", "35b"),
    spec("hla", "HLA typing (T1K)", "08"),
    spec("prs", "Polygenic risk scores", "25"),
    spec(
        "ancestry_fine",
        "Fine-grained ancestry (ADMIXTURE, RFMix)",
        "26b",
    ),
    spec("manta", "Structural variants (Manta)", "04"),
    spec("delly", "Structural variants (Delly)", "19"),
    spec("cnvpytor", "Copy number (CNVpytor)", "18"),
    spec("sv_consensus", "SV consensus merge", "22"),
    spec(
        "annotsv",
        "Structural variants, class 4 and 5 (AnnotSV)",
        "05",
    ),
    spec("expansions", "Repeat expansions (ExpansionHunter)", "09"),
    spec("telomere", "Telomere content (TelomereHunter)", "10"),
    spec("roh", "Runs of homozygosity", "11"),
    spec("haplogroup", "Mitochondrial haplogroup", "12"),
    spec("y_haplogroup", "Y-chromosome haplogroup (Yleaf)", "37"),
    spec("y_tree", "Y placement on the current public tree", "37"),
    spec("y_str", "Y-STR markers (HipSTR)", "45"),
    spec(
        "consumer_export",
        "Raw-data export (23andMe v5 layout)",
        "40",
    ),
    spec("bloodtype", "Blood groups (predicted)", "41"),
    spec("lpa", "Lp(a): LPA KIV-2 repeat count", "46"),
    spec("vntyper", "MUC1 VNTR frameshifts (VNtyper)", "52"),
    spec("locityper", "Locityper second opinion at hard loci", "53"),
    spec("mito", "Mitochondrial variants", "20"),
    spec("traits", "Traits panel", "42"),
    spec("archaic", "Archaic ancestry (IBDmix)", "43"),
    spec("ancient", "Ancient populations (AADR PCA)", "44"),
    spec("viral", "Viral DNA in the unmapped reads", "47"),
    spec("exomiser", "Symptom-driven ranking (Exomiser)", "48"),
    spec("chip", "Clonal haematopoiesis (CHIP)", "49"),
    spec("mei", "Mobile element insertions (MEGAnE)", "50"),
    spec("mosaic", "Mosaic chromosomal alterations (MoChA)", "51"),
];

const fn spec(key: &'static str, title: &'static str, step: &'static str) -> Spec {
    Spec { key, title, step }
}

/// A file a section found: its path and its path relative to the folder.
#[derive(Debug, Clone)]
pub struct Found {
    pub path: String,
    pub rel: String,
}

/// What a reader sees: the sample, its folder, the assets and one walk of
/// the folder.
pub struct Ctx<'a> {
    pub sample: &'a str,
    pub dir: &'a str,
    pub assets: &'a str,
    /// The folder of `hla_alleles.tsv` and `str_intermediate_ranges.tsv`.
    pub tables: &'a str,
    pub index: &'a FolderIndex,
}

impl Ctx<'_> {
    /// The sample name as a file pattern matches it literally.
    pub fn s(&self) -> String {
        pycompat::mtime::escape(self.sample)
    }

    /// `first_existing(d, patterns)`.
    pub fn find(&self, pats: &[String]) -> Option<Found> {
        self.index.first(pats).map(|rel| Found {
            path: crate::join(self.dir, rel),
            rel: rel.to_owned(),
        })
    }

    /// `first_existing` of one pattern.
    pub fn find1(&self, pat: String) -> Option<Found> {
        self.find(&[pat])
    }

    /// `os.path.join(d, rel)`.
    pub fn path(&self, rel: &str) -> String {
        crate::join(self.dir, rel)
    }

    /// `os.path.join(os.path.dirname(found), name)` when it is a file.
    pub fn beside(&self, f: &Found, name: &str) -> Option<Found> {
        let dir = |p: &str| p.rsplit_once('/').map_or(String::new(), |x| x.0.to_owned());
        let path = crate::join(&dir(&f.path), name);
        crate::is_file(&path).then(|| Found {
            rel: crate::join(&dir(&f.rel), name),
            path,
        })
    }

    /// A table of the `bin/` folder ([`util::default_tables`]).
    pub fn table(&self, name: &str) -> String {
        crate::join(self.tables, name)
    }

    /// `assets/<name>`.
    pub fn asset(&self, name: &str) -> String {
        crate::join(self.assets, name)
    }
}

/// What a section reader gives: the file, its data, and for the sections
/// the post pass reads, what it needs, typed.
pub type Out = Result<Option<(Found, Value)>, String>;

/// The states and data of group b the post pass reads. Empty under
/// `--only c`, which reads no group b section.
#[derive(Debug, Default, Clone)]
pub struct Others {
    pub cpic: Option<State>,
    pub traits: Option<State>,
    /// pypgx's CYP2D6 call and its section state, when that section is ok
    /// or stale.
    pub pypgx: Option<(State, String)>,
    /// Cyrius's genotype and Filter (either may be absent), when that
    /// section is ok or stale.
    pub cyrius: Option<(State, Option<String>, Option<String>)>,
}

/// Reads every section of group c, in [`GROUP_C`] order.
pub fn read_all(c: &Ctx<'_>) -> Vec<Out> {
    let readers: [fn(&Ctx<'_>) -> Out; 32] = [
        tables::couple,
        cn::parascopy,
        cn::paralogs,
        hla::hla,
        tables::prs,
        tables::ancestry_fine,
        sv::manta,
        sv::delly,
        sv::cnvpytor,
        sv::sv_consensus,
        sv::annotsv,
        repeats::expansions,
        tables::telomere,
        roh::section,
        tables::haplogroup,
        tables::y_haplogroup,
        tables::y_tree,
        tables::y_str,
        tables::consumer_export,
        tables::bloodtype,
        tables::lpa,
        tables::vntyper,
        tables::locityper,
        chrm::mito,
        traits::traits,
        tables::archaic,
        tables::ancient,
        tables::viral,
        tables::exomiser,
        tables::chip,
        tables::mei,
        tables::mosaic,
    ];
    // The readers are independent: one thread each, over a small pool, so a
    // big input (a VCF, a table of thousands of rows) never holds up the
    // others. A panic stays in its section.
    let jobs = std::thread::available_parallelism().map_or(4, |n| n.get().min(8));
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<Out>>> = (0..readers.len())
        .map(|_| std::sync::Mutex::new(None))
        .collect();
    std::thread::scope(|scope| {
        for _ in 0..jobs {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(r) = readers.get(i) else { break };
                    let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| r(c)))
                        .unwrap_or_else(|_| Err("internal error: the reader stopped".to_owned()));
                    if let Ok(mut slot) = slots[i].lock() {
                        *slot = Some(out);
                    }
                }
            });
        }
    });
    slots
        .into_iter()
        .map(|m| {
            m.into_inner()
                .ok()
                .flatten()
                .unwrap_or_else(|| Err("internal error: the reader stopped".to_owned()))
        })
        .collect()
}

/// The post pass over the read sections (`sections` in [`GROUP_C`]
/// order, each with its state and data): the repeat readings with the
/// traits and CPIC states, each Locityper locus's primary call, and the
/// chrY filter of CNVpytor's count. Each part runs on its own: one that
/// returns an error is returned as its section's key and why, and the
/// others still run. The parts read with defaults, so a row of an
/// unexpected shape is no error. A panic is not caught: the release binary
/// aborts on one (`panic = "abort"`).
pub fn post(
    c: &Ctx<'_>,
    sections: &mut [(Spec, State, Option<&mut Data>)],
    others: &Others,
    declared_sex: Option<&str>,
) -> Vec<(&'static str, String)> {
    let state = |key: &str, secs: &[(Spec, State, Option<&mut Data>)]| {
        secs.iter()
            .find(|(s, ..)| s.key == key)
            .map_or(State::Missing, |x| x.1)
    };
    let traits_state = others.traits.unwrap_or_else(|| state("traits", sections));
    let cpic_state = others.cpic.unwrap_or(State::Ok);
    let hla_state = state("hla", sections);
    let hla_data = sections
        .iter()
        .find(|(s, ..)| s.key == "hla")
        .and_then(|x| x.2.as_deref())
        .and_then(json)
        .cloned();
    let mut failed = Vec::new();
    for (spec, st, data) in sections.iter_mut() {
        let live = matches!(st, State::Ok | State::Stale);
        let Some(Data::Json(v)) = data.as_deref_mut() else {
            continue;
        };
        let part = || -> Result<(), String> {
            match spec.key {
                "expansions" if live => repeats::post(c, v, traits_state, cpic_state)
                    .map_err(|e| format!("the repeat readings could not be made: {e}")),
                "locityper" if live => {
                    tables::locityper_post(v, hla_state, hla_data.as_ref(), others);
                    Ok(())
                }
                "cnvpytor" => {
                    sv::cnvpytor_sex_filter(v, declared_sex);
                    Ok(())
                }
                _ => Ok(()),
            }
        };
        failed.extend(contain(spec.key, part));
    }
    failed
}

/// Runs one part of the post pass: `None` when it went through, else its
/// section's key and the error it returned.
pub fn contain(
    key: &'static str,
    part: impl FnOnce() -> Result<(), String>,
) -> Option<(&'static str, String)> {
    part().err().map(|e| (key, e))
}

fn json(d: &Data) -> Option<&Value> {
    match d {
        Data::Json(v) => Some(v),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
