//! Both reports point to where a haplogroup can be explored further (steps
//! 12 and 37). It replaces `tests/test_haplogroup_explore.py`.
//!
//! Each case writes the haplogroup files of one made-up sample S into an
//! output folder, collects the summary and renders both reports:
//!
//! 1. a male sample (Yleaf table with R-M269, haplogrep3 with H5a7): the Y
//!    and the mtDNA links are built from those names, in the text and the
//!    HTML report, with the note that nothing is uploaded;
//! 2. Yleaf's paragroup form R-M269*(xR-L23): the links use R-M269;
//! 3. a female sample (no Yleaf table): no Y links; the mtDNA links stay;
//! 4. no call (Hg NA, haplogroup ''): no links at all;
//! 5. a branch name with characters a URL cannot carry raw is encoded;
//! 6. an address that would close the link's attribute and open a tag is
//!    escaped in the HTML report, so the page never runs it.
//!
//! The branch names are the fixture's and public examples, not anyone's
//! result.

use std::path::{Path, PathBuf};

use og_render::model::Summary;
use og_render::{html, text};
use pycompat::json::{Dump, dumps};

const Y_HEAD: &str =
    "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn put(path: &Path, text: &str) {
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    std::fs::write(path, text).unwrap();
}

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

/// A sample folder with haplogrep3's call and, for a male, Yleaf's table.
fn sample(name: &str, y_hg: Option<&str>, mt_hg: &str) -> Tmp {
    let d = std::env::temp_dir().join(format!("og-ancestry-hg-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    put(
        &d.join("mito/S_haplogroup.txt"),
        &format!("\"SampleID\"\t\"Haplogroup\"\t\"Rank\"\n\"S\"\t\"{mt_hg}\"\t\"0.98\"\n"),
    );
    if let Some(hg) = y_hg {
        put(
            &d.join("y_haplogroup/S_y_haplogroup.txt"),
            &format!("{Y_HEAD}S_sorted\t{hg}\tM269\t100\t42\t0.97\t1\t1\t0.97\n"),
        );
    }
    Tmp(d)
}

/// The summary og summary collect writes for S, as text.
fn collect(d: &Path) -> String {
    let input = og_summary::collect::Input {
        sample: "S".to_owned(),
        sample_dir: d.display().to_string(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: og_summary::collect::Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    dumps(
        &og_summary::collect::post::collect(&input).unwrap(),
        &Dump::default(),
    )
    .unwrap()
}

/// Both reports of a summary.
fn render(summary: &str) -> (String, String) {
    let s = Summary::parse(summary).unwrap();
    (text::report(&s).unwrap(), html::report(&s).unwrap())
}

fn reports(d: &Path) -> (String, String) {
    render(&collect(d))
}

#[test]
fn a_male_sample_gets_the_y_and_the_mtdna_links() {
    let (txt, page) = reports(&sample("male", Some("R-M269"), "H5a7"));
    for url in [
        "https://discover.familytreedna.com/y-dna/R-M269",
        "https://www.yfull.com/tree/R-M269/",
        "https://discover.familytreedna.com/mtdna/H5a7",
        "https://www.yfull.com/mtree/H5a7/",
    ] {
        assert!(txt.contains(url), "{url} in the text report:\n{txt}");
        assert!(
            page.contains(&format!("<a href=\"{url}\">")),
            "{url} as a link in the HTML report"
        );
    }
    assert!(txt.contains("docs/37-y-haplogroup.md#explore-your-branch-further"));
    assert!(txt.contains("docs/12-mito-haplogroup.md#explore-your-branch-further"));
    assert!(txt.contains("Nothing is uploaded by this pipeline"));
    assert!(page.contains("Nothing is uploaded by this pipeline"));
}

#[test]
fn a_paragroup_links_its_named_branch() {
    let (txt, _) = reports(&sample("para", Some("R-M269*(xR-L23)"), "H5a7"));
    assert!(txt.contains("https://www.yfull.com/tree/R-M269/"), "{txt}");
    assert!(!txt.contains("y-dna/R-M269*"), "{txt}");
    assert!(!txt.contains("%2A"), "{txt}");
}

#[test]
fn a_female_sample_has_no_y_link_and_keeps_the_mtdna_ones() {
    let (txt, page) = reports(&sample("female", None, "H5a7"));
    for t in [&txt, &page] {
        assert!(!t.contains("y-dna/"), "{t}");
        assert!(!t.contains("yfull.com/tree/"), "{t}");
    }
    assert!(txt.contains("https://discover.familytreedna.com/mtdna/H5a7"));
    assert!(page.contains("https://www.yfull.com/mtree/H5a7/"));
}

#[test]
fn no_call_has_no_link() {
    let (txt, page) = reports(&sample("nocall", Some("NA"), ""));
    assert!(!txt.contains("discover.familytreedna.com"), "{txt}");
    assert!(!txt.contains("yfull.com"), "{txt}");
    assert!(!page.contains("discover.familytreedna.com"));
    assert!(
        txt.contains("insufficient markers"),
        "the Y section still says why:\n{txt}"
    );
}

#[test]
fn a_plus_in_the_branch_is_encoded() {
    let (txt, _) = reports(&sample("esc", None, "L3e2b+152"));
    assert!(txt.contains("mtdna/L3e2b%2B152"), "{txt}");
    assert!(txt.contains("mtree/L3e2b%2B152/"), "{txt}");
}

#[test]
fn a_planted_address_cannot_open_a_tag_in_the_html_report() {
    let summary = collect(&sample("planted", Some("R-M269"), "H5a7"));
    let planted = r#"https://x.example/\"><script>alert(1)</script><a href=\""#;
    let mut s = summary.clone();
    for url in [
        "https://www.yfull.com/mtree/H5a7/",
        "https://www.yfull.com/tree/R-M269/",
    ] {
        assert!(s.contains(url), "the summary links {url}");
        s = s.replace(url, planted);
    }
    let (txt, page) = render(&s);
    assert!(
        !page.contains("<script"),
        "a planted address opened a tag in the HTML report"
    );
    assert!(!page.contains("x.example/\"><"), "it ended the attribute");
    assert_eq!(
        page.matches("&quot;&gt;&lt;script&gt;alert(1)&lt;/script&gt;")
            .count(),
        4,
        "escaped in the address and the text of both links"
    );
    assert!(
        txt.contains("https://x.example/\"><script>"),
        "the text report shows it as written"
    );
}

#[test]
fn an_address_that_is_not_a_web_address_is_shown_as_text() {
    let summary = collect(&sample("scheme", Some("R-M269"), "H5a7"));
    let url = "https://www.yfull.com/mtree/H5a7/";
    assert!(summary.contains(url), "the summary links {url}");
    let s = summary.replace(url, "javascript:alert(1)");
    let (txt, page) = render(&s);
    assert!(
        !page.contains("href=\"javascript:"),
        "a javascript: address became a link in the HTML report"
    );
    assert!(
        page.contains("<span class=\"value\">javascript:alert(1)</span>"),
        "it is shown as text"
    );
    assert!(
        page.contains("<a href=\"https://www.yfull.com/tree/R-M269/\">"),
        "a web address stays a link"
    );
    assert!(txt.contains("javascript:alert(1)"), "{txt}");
}
