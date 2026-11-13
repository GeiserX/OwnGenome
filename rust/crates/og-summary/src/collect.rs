//! `og summary collect`: one sample folder's step outputs read into
//! summary.json, the file both reports, the live site and the run delta
//! read.
//!
//! The port comes in groups (`--only`), each checked against the Python
//! collector (`bin/collect_summary.py`) section by section before anything
//! switches to it. Group `a` is the model and the first seven sections:
//! coverage, callability, sex_check, sample_qc, mosaic_sex_loss,
//! read_count and variants, with the run status, the manifest and each
//! section's state. Group `c` is copy number, HLA, the polygenic scores,
//! ancestry, structural variants, repeats, telomeres, runs of
//! homozygosity, the haplogroups, the opt-in step tables, chrM, traits with
//! the APOL1 G2 deletion, and the couple screen ([`group_c`]).
//!
//! How it reads a folder: one walk of the folder ([`scan::FolderIndex`])
//! instead of a `glob` per candidate path; each input read once (the
//! step 33 table feeds two sections); the sample VCF, the one big input,
//! read on its own thread with BGZF inflated on many cores while the other
//! sections are read. A reader's error stays in its section, which is
//! then `unreadable` with a note that names the file.

pub mod group_c;
pub mod model;
pub mod scan;
pub mod sections;
pub mod time;

use crate::args::{self, Kind, Opt};
use model::{GROUP_A, Manifest, RunStatus, Section, Spec, Summary};
use pycompat::json::{Dump, dumps};
use scan::FolderIndex;
use sections::{Callability, Coverage, Data, Mosaic, ReadCount, SampleQc, SexCheck, Variants};
use std::io::Write;
use std::path::Path;

const USAGE: &str = "usage: og summary collect --sample SAMPLE --sample-dir SAMPLE_DIR --only {a,c}
                          --out OUT [--declared-sex DECLARED_SEX] [--assets ASSETS]
";

const OPTS: &[Opt] = &[
    Opt("--sample", Kind::Value, true),
    Opt("--sample-dir", Kind::Value, true),
    Opt("--declared-sex", Kind::Value, false),
    Opt("--out", Kind::Value, true),
    Opt("--assets", Kind::Value, false),
    Opt("--only", Kind::Choice(&["a", "c"]), true),
];

/// What one collect reads.
#[derive(Debug, Clone)]
pub struct Input {
    pub sample: String,
    pub sample_dir: String,
    /// The samplesheet's sex, for a folder whose run status and manifest
    /// record none.
    pub declared_sex: Option<String>,
    /// The repository's `assets/` folder.
    pub assets: String,
    /// The group of sections to read: `a` or `c`.
    pub group: Group,
    /// The folder of the two tables group c reads from `bin/`.
    pub tables: String,
}

/// The group of sections a collect reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Group {
    A,
    C,
}

/// `og summary collect`: exit 0 when summary.json is written (whatever
/// its sections say), 1 when it cannot be written, 2 for a command line it
/// cannot read.
pub fn cli(argv: &[String], out: &mut impl Write, err: &mut impl Write) -> u8 {
    let a = match args::parse(argv, OPTS) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og summary collect: error: {e}");
            return 2;
        }
    };
    let input = Input {
        sample: a.text("--sample").to_owned(),
        sample_dir: a.text("--sample-dir").to_owned(),
        declared_sex: a.get("--declared-sex").map(str::to_owned),
        assets: a.get("--assets").map_or_else(default_assets, str::to_owned),
        tables: group_c::util::default_tables(),
        group: if a.text("--only") == "c" {
            Group::C
        } else {
            Group::A
        },
    };
    // A part of the sections must never replace the sample's complete
    // summary.json, so there is no default path while `--only` exists.
    let target = a.text("--out").to_owned();
    let summary = match collect(&input) {
        Ok(s) => s,
        Err(e) => return crate::finish(Err(e), err),
    };
    let body = match dumps(
        &summary.to_json(),
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    ) {
        Ok(b) => b + "\n",
        Err(e) => return crate::finish(Err(e.to_string()), err),
    };
    let r = crate::write_file(&target, body.as_bytes());
    if r.is_ok() {
        let _ = writeln!(out, "Summary written: {target}");
    }
    crate::finish(r, err)
}

/// The `assets/` folder of the checkout the binary sits in (`bin/.og/og-*`
/// and `rust/target/*/og` both do), else the one beside the binary's
/// folder, as the Python found it beside `bin/`.
fn default_assets() -> String {
    let exe = std::env::current_exe().unwrap_or_default();
    exe.ancestors()
        .skip(1)
        .map(|d| d.join("assets"))
        .find(|p| p.join("callability_genes.tsv").is_file())
        .unwrap_or_else(|| {
            exe.parent()
                .and_then(Path::parent)
                .unwrap_or(Path::new("."))
                .join("assets")
        })
        .display()
        .to_string()
}

/// Reads the folder into the summary of group a. An error is a run status
/// or a manifest that is there and cannot be read: without it a result of
/// an earlier run would be called current, so no summary is written.
pub fn collect(input: &Input) -> Result<Summary, String> {
    let dir = Path::new(&input.sample_dir);
    let index = FolderIndex::walk(dir);
    let run =
        read_optional(&dir.join("logs").join("run_status.tsv"))?.map(|t| RunStatus::parse(&t));
    let manifest = read_optional(&dir.join("run_manifest.tsv"))?.map(|t| Manifest::parse(&t));
    let s = pycompat::mtime::escape(&input.sample);
    let find = |pats: &[String]| index.first(pats).map(str::to_owned);
    let file = |rel: String| -> (String, String) { (crate::join(&input.sample_dir, &rel), rel) };
    if input.group == Group::C {
        // The declared sex the chrY filter of CNVpytor reads comes from
        // step 33's table when the run status and manifest have none.
        let qc = find(&[format!("qc/{s}_sample_qc.tsv")]).map(file);
        let qc_declared = qc.and_then(|(p, _)| {
            let rows = sections::read_tsv(&p).ok()?;
            SampleQc::read(&p, &rows)
                .ok()?
                .declared_sex()
                .map(str::to_owned)
        });
        let declared_sex = declared(
            input,
            run.as_ref(),
            manifest.as_ref(),
            qc_declared.as_deref(),
        );
        let sections = read_group_c(input, &index, run.as_ref(), declared_sex.as_deref());
        return Ok(Summary {
            sample: input.sample.clone(),
            generated_utc: time::utc_minute(now()).unwrap_or_default(),
            run,
            manifest,
            declared_sex,
            sections,
        });
    }

    let variants_at = find(&[format!("vcf/{s}.vcf.gz")]);
    let (variants, rest) = std::thread::scope(|scope| {
        let v = scope.spawn(|| {
            read_at(variants_at.map(file), |p| {
                Variants::read(p).map(Data::Variants)
            })
        });
        let rest = read_rest(input, &find, &file);
        (
            v.join()
                .unwrap_or_else(|_| Err("internal error: the variant reader stopped".into())),
            rest,
        )
    });
    let [
        coverage,
        callability,
        sex_check,
        sample_qc,
        mosaic,
        read_count,
    ] = rest;
    let reads = [
        coverage,
        callability,
        sex_check,
        sample_qc,
        mosaic,
        read_count,
        variants,
    ];
    let sections: Vec<Section> = GROUP_A
        .iter()
        .zip(reads)
        .map(|(spec, r): (&Spec, _)| Section::new(*spec, r, run.as_ref()))
        .collect();

    let qc_declared = sections
        .iter()
        .find(|s| s.spec.key == "sample_qc")
        .and_then(|s| match &s.data {
            Some(Data::SampleQc(q)) => q.declared_sex(),
            _ => None,
        });
    let declared_sex = declared(input, run.as_ref(), manifest.as_ref(), qc_declared);
    Ok(Summary {
        sample: input.sample.clone(),
        generated_utc: time::utc_minute(now()).unwrap_or_default(),
        run,
        manifest,
        declared_sex,
        sections,
    })
}

fn now() -> f64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0.0, |d| d.as_secs_f64())
}

/// The declared sex: the command line's, the run status's, the
/// manifest's, then step 33's; the first that is not empty.
fn declared(
    input: &Input,
    run: Option<&RunStatus>,
    manifest: Option<&Manifest>,
    qc: Option<&str>,
) -> Option<String> {
    [
        input.declared_sex.as_deref(),
        run.and_then(|r| r.declared_sex.as_deref()),
        manifest.and_then(|m| m.value("run", "declared_sex")),
        qc,
    ]
    .into_iter()
    .flatten()
    .find(|v| !v.is_empty())
    .map(str::to_owned)
}

/// Reads the sections of group c, then runs the parts of the post pass
/// that change them.
fn read_group_c(
    input: &Input,
    index: &FolderIndex,
    run: Option<&RunStatus>,
    declared_sex: Option<&str>,
) -> Vec<Section> {
    let c = group_c::Ctx {
        sample: &input.sample,
        dir: &input.sample_dir,
        assets: &input.assets,
        tables: &input.tables,
        index,
    };
    let reads = group_c::read_all(&c);
    let mut sections: Vec<Section> = group_c::GROUP_C
        .iter()
        .zip(reads)
        .map(|(spec, r)| {
            let r = r.and_then(|found| match found {
                None => Ok(None),
                Some((f, v)) => with_mtime(&f.path, f.rel, Data::Json(v)),
            });
            Section::new(*spec, r, run)
        })
        .collect();
    let mut view: Vec<(Spec, model::State, Option<&mut Data>)> = sections
        .iter_mut()
        .map(|s| (s.spec, s.state, s.data.as_mut()))
        .collect();
    group_c::post(&c, &mut view, &group_c::Others::default(), declared_sex);
    sections
}

/// What one reader gave: the file (relative), its mtime and its data.
type Read = Result<Option<(String, f64, Data)>, String>;

/// Reads the six sections other than variants, in order.
fn read_rest(
    input: &Input,
    find: &dyn Fn(&[String]) -> Option<String>,
    file: &dyn Fn(String) -> (String, String),
) -> [Read; 6] {
    let s = pycompat::mtime::escape(&input.sample);
    let coverage = read_at(
        find(&[
            format!("mosdepth/{s}.mosdepth.summary.txt"),
            format!("coverage/{s}.mosdepth.summary.txt"),
        ])
        .map(file),
        |p| Coverage::read(p).map(Data::Coverage),
    );
    let callability = read_at(
        find(&[
            format!("mosdepth/{s}_callability.bed.gz"),
            format!("coverage/{s}_callability.bed.gz"),
        ])
        .map(file),
        |p| {
            let catalog = sections::read_catalog(&input.assets)?;
            Callability::read(&catalog, p).map(Data::Callability)
        },
    );
    let sex_check = read_at(
        find(&[
            "indexcov/indexcov-indexcov.ped".to_owned(),
            "indexcov/*-indexcov.ped".to_owned(),
        ])
        .map(file),
        |p| SexCheck::read(p).map(Data::SexCheck),
    );
    // Step 33's table feeds two sections; it is read once.
    let qc_at = find(&[format!("qc/{s}_sample_qc.tsv")]).map(file);
    let qc_rows = qc_at.as_ref().map(|(p, _)| sections::read_tsv(p));
    let rows = |p: &str| qc_rows.clone().unwrap_or_else(|| sections::read_tsv(p));
    let sample_qc = read_at(qc_at.clone(), |p| {
        SampleQc::read(p, &rows(p)?).map(Data::SampleQc)
    });
    let mosaic = match qc_at {
        None => Ok(None),
        Some((p, rel)) => match rows(&p).and_then(|r| Mosaic::read(&p, &r)) {
            Err(e) => Err(e),
            Ok(None) => Ok(None),
            Ok(Some(m)) => with_mtime(&p, rel, Data::MosaicSexLoss(m)),
        },
    };
    let read_count = read_at(find(&[format!("qc/{s}_read_count.tsv")]).map(file), |p| {
        ReadCount::read(p).map(Data::ReadCount)
    });
    [
        coverage,
        callability,
        sex_check,
        sample_qc,
        mosaic,
        read_count,
    ]
}

/// Runs a reader on the file found, if one was.
fn read_at(
    found: Option<(String, String)>,
    reader: impl FnOnce(&str) -> Result<Data, String>,
) -> Read {
    let Some((path, rel)) = found else {
        return Ok(None);
    };
    let data = reader(&path)?;
    with_mtime(&path, rel, data)
}

fn with_mtime(path: &str, rel: String, data: Data) -> Read {
    let t = pycompat::mtime::mtime(path).map_err(|e| format!("{}: {e}", crate::basename(path)))?;
    Ok(Some((rel, t, data)))
}

/// The text of a file that may not be there: `None` when it is no file,
/// an error when it is one and cannot be read.
fn read_optional(p: &Path) -> Result<Option<String>, String> {
    let Some(p) = p.to_str() else {
        return Ok(None);
    };
    if !crate::is_file(p) {
        return Ok(None);
    }
    crate::read_text(p).map(Some)
}

#[cfg(test)]
mod tests;
