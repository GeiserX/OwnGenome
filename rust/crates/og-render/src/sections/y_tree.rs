//! The second Y placement (step 37, `og y tree-place`): the branch on a
//! public tree release newer than Yleaf's, the SNP counts on the path to
//! it, and one row per level of the walk.

use crate::model::{Card, Cell, Inline, Layout, Line, Malformed, State, Summary, Table};
use crate::notes::{EXPLORE_NOTE, Y_TREE_NOTE};
use crate::sections::ancestry::{
    Align, explore_html, get_text, list, need_align, need_text, py_or, shown,
};
use crate::{html, text};

const PATH: &str = "sections.y_tree.data";
const ROW: &str = "sections.y_tree.data.path[]";

/// The two lists of ancestral calls a reader should see, with their labels.
const ANCESTRAL: [(&str, &str); 2] = [
    (
        "Ancestral inside the terminal branch",
        "ancestral_in_terminal",
    ),
    (
        "Ancestral on the path, outvoted by the SNPs below it",
        "outvoted_ancestral",
    ),
];

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "y_tree", Some("Y Placement on the Current Public Tree"))?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = s.section("y_tree").data;
    let n = |k: &str| need_text(d, k, PATH);
    w.push(format!(
        "  Branch: {} (level {}, {})",
        n("haplogroup")?,
        n("level")?,
        n("tree")?
    ));
    w.push(format!(
        "  Path SNPs: {} derived, {} ancestral, {} no call; {} SNPs called at depth {} or more",
        n("path_derived")?,
        n("path_ancestral")?,
        n("path_no_call")?,
        n("snps_called")?,
        n("min_depth")?
    ));
    for (label, key) in ANCESTRAL {
        if shown(d, key, PATH)?.is_some() {
            w.push(format!("  {label}: {}", n(key)?));
        }
    }
    for x in list(d, "path", PATH)? {
        let a = |k: &str, width: usize, side: Align| need_align(&x, k, width, side, ROW);
        w.push(format!(
            "  {}  {} {}D {}A {}N",
            a("Level", 3, Align::Right)?,
            a("Branch", 22, Align::Left)?,
            a("Derived", 3, Align::Right)?,
            a("Ancestral", 3, Align::Right)?,
            a("No_call", 4, Align::Right)?
        ));
    }
    w.push(format!("  {Y_TREE_NOTE}"));
    w.push(String::new());
    Ok(w)
}

/// The HTML report's card, when step 37 placed the sample on the tree.
pub fn card(s: &Summary) -> Result<Option<Card>, Malformed> {
    let sec = s.section("y_tree");
    if sec.state == State::Missing {
        return Ok(None);
    }
    let mut c = Card {
        title: "Y Placement on the Current Public Tree".to_owned(),
        full: true,
        ..Card::default()
    };
    if sec.state.ran() {
        let d = sec.data;
        let g = |k: &str| get_text(d, k, ".", PATH);
        for (label, key) in [
            ("Branch", "haplogroup"),
            ("Level", "level"),
            ("Tree", "tree"),
        ] {
            c.body
                .push(Line::Stat(label.to_owned(), Inline::Text(g(key)?)));
        }
        c.body.push(Line::Para(format!(
            "Path SNPs: {} derived, {} ancestral, {} no call; {} SNPs called at depth {} or more",
            g("path_derived")?,
            g("path_ancestral")?,
            g("path_no_call")?,
            g("snps_called")?,
            g("min_depth")?
        )));
        for (label, key) in ANCESTRAL {
            if let Some(v) = shown(d, key, PATH)? {
                c.body
                    .push(Line::Para(format!("{label}: {}", py_or(Some(v), "")?)));
            }
        }
        let mut rows = Vec::new();
        for x in list(d, "path", PATH)? {
            rows.push(
                [
                    "Level",
                    "Branch",
                    "Derived",
                    "Ancestral",
                    "No_call",
                    "Ancestral_SNPs",
                ]
                .iter()
                .map(|k| need_text(&x, k, ROW).map(Cell::Text))
                .collect::<Result<Vec<_>, Malformed>>()?,
            );
        }
        if !rows.is_empty() {
            c.body.push(Line::Table(Table {
                head: [
                    "Level",
                    "Branch",
                    "Derived",
                    "Ancestral",
                    "No call",
                    "Ancestral SNPs",
                ]
                .map(str::to_owned)
                .to_vec(),
                rows,
                layout: Layout::OneLine,
                ..Table::default()
            }));
        }
        let (links, note) = explore_html(d, PATH)?;
        c.body.extend(links);
        if note {
            c.caveats.push(EXPLORE_NOTE);
        }
        c.caveats.push(Y_TREE_NOTE);
    }
    Ok(Some(html::card(&sec, c)?))
}
