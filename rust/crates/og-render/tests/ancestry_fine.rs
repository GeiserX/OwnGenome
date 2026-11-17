//! Both reports of step 26b's table (`og ancestry-fine report` writes it):
//! the groups with their bootstrap range, the group the model leaves out,
//! the chromosome RFMix failed on, the painted shares and the painting.
//! It replaces the render part of `tests/test_ancestry_fine_report.py`;
//! og-ancestry's own tests check the table the report command writes.

use std::path::{Path, PathBuf};

use og_render::model::Summary;
use og_render::{html, text};
use pycompat::json::{Dump, dumps};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// The table `og ancestry-fine report` wrote for the Python test's inputs.
const TABLE: &str = "level\tcode\tgroup\tfraction\tse\tlow\thigh
continent\tWEA\tWest Eurasia (Europe, Middle East, North Africa)\t0.9000\t0.0332\t0.8025\t0.8975
continent\tEAS\tEast Asia\t0.0500\t0.0000\t0.0500\t0.0500
continent\tAFR\tAfrica\t0.0500\t0.0332\t0.0525\t0.1475
region\tIBE\tIberia\t0.6000\t.\t.\t.
region\tITA\tItaly\t0.4000\t.\t.\t.
painting\tWEA\tWest Eurasia (Europe, Middle East, North Africa)\t0.7000\t.\t.\t.
painting\tAFR\tAfrica\t0.3000\t.\t.\t.
left_out\tOCE\tOceania (continental group)\t.\t.\t.\t.
left_out\tOCE\tOceania (regional group)\t.\t.\t.\t.
not_painted\tchr19\tRFMix failed\t.\t.\t.\t.
";

/// A scratch folder that is removed when the test ends, also when it fails.
struct Tmp(PathBuf);

impl std::ops::Deref for Tmp {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for Tmp {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<std::ffi::OsStr> for Tmp {
    fn as_ref(&self) -> &std::ffi::OsStr {
        self.0.as_os_str()
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn sample(tag: &str, painting: bool) -> Tmp {
    let d = std::env::temp_dir().join(format!("og-ancestry-fine-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("ancestry_fine")).unwrap();
    std::fs::write(d.join("ancestry_fine/S_ancestry_fine.tsv"), TABLE).unwrap();
    if painting {
        std::fs::write(
            d.join("ancestry_fine/S_painting.png"),
            b"\x89PNG\r\n\x1a\nnot a real picture",
        )
        .unwrap();
    }
    Tmp(d)
}

fn reports(d: &Path) -> (String, String) {
    let input = og_summary::collect::Input {
        sample: "S".to_owned(),
        sample_dir: d.display().to_string(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: og_summary::collect::Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    let summary = dumps(
        &og_summary::collect::post::collect(&input).unwrap(),
        &Dump::default(),
    )
    .unwrap();
    let s = Summary::parse(&summary).unwrap();
    (text::report(&s).unwrap(), html::report(&s).unwrap())
}

#[test]
fn both_reports_name_the_group_left_out_and_the_chromosome_not_painted() {
    let (txt, page) = reports(&sample("full", true));
    for doc in [&txt, &page] {
        assert!(
            doc.contains(
                "Left out of the model: Oceania (continental and regional group). The panel has \
                 fewer than 10"
            ),
            "{doc}"
        );
        assert!(doc.contains("not painted: RFMix failed on chr19"), "{doc}");
    }
    assert!(txt.contains("Fine-Grained Ancestry (ADMIXTURE, RFMix)"));
    assert!(
        txt.contains("    West Eurasia (Europe, Middle East, North Africa)     90% (80% to 90%)"),
        "{txt}"
    );
    assert!(
        txt.contains(
            "  Chromosome painting, share of the painted genome: West Eurasia (Europe, Middle \
             East, North Africa) 70%, Africa 30%; not painted: RFMix failed on chr19"
        ),
        "{txt}"
    );
    assert!(page.contains("<h2>Fine-Grained Ancestry</h2>"));
    assert!(page.contains(
        "<tr><td>West Eurasia (Europe, Middle East, North Africa)</td><td>90%</td><td>80% to \
         90%</td></tr>"
    ));
    assert!(
        page.contains("<tr><td>Iberia</td><td>60%</td><td>.</td></tr>"),
        "a level without bootstrap replicates has no range"
    );
    assert!(
        page.contains("src=\"data:image/png;base64,iVBORw0KGg"),
        "the painting"
    );
    assert!(
        page.contains("&#9632;</span> not painted</p>"),
        "its legend"
    );
}

#[test]
fn without_a_picture_the_card_has_no_image() {
    let (_, page) = reports(&sample("nopng", false));
    assert!(page.contains("<h2>Fine-Grained Ancestry</h2>"));
    assert!(!page.contains("data:image/png;base64,"));
    assert!(
        !page.contains("not painted</p>"),
        "no legend without the painting"
    );
}
