//! The report's chapters: each one opens with the question it answers and
//! holds sections, and each section shows one or more summary keys. Both
//! reports and the live site lay out the results by this one table.

use og_summary::collect::{group_b::GROUP_B, group_c::GROUP_C, model::GROUP_A, model::Spec};

use crate::model::{Colour, HtmlPart, Malformed, State, Summary};
use crate::sections;

/// One section of a chapter: its anchor, its title in the contents and the
/// summary keys it shows.
#[derive(Debug, Clone, Copy)]
pub struct Sec {
    pub id: &'static str,
    pub title: &'static str,
    pub keys: &'static [&'static str],
}

/// One chapter: its anchor (`ch-<id>`), title, question and sections.
#[derive(Debug, Clone, Copy)]
pub struct Chapter {
    pub id: &'static str,
    pub title: &'static str,
    pub question: &'static str,
    pub sections: &'static [Sec],
}

const fn sec(id: &'static str, title: &'static str, keys: &'static [&'static str]) -> Sec {
    Sec { id, title, keys }
}

/// The chapters, in the order both reports show them.
pub const CHAPTERS: [Chapter; 7] = [
    Chapter {
        id: "quality",
        title: "Quality",
        question: "Can I trust this data?",
        sections: &[
            sec(
                "qc",
                "Quality control",
                &["coverage", "sex_check", "sample_qc", "read_count"],
            ),
            sec(
                "callability",
                "Genes this genome could not read reliably",
                &["callability"],
            ),
            sec("variants", "Variant calling", &["variants"]),
        ],
    },
    Chapter {
        id: "clinical",
        title: "Clinical",
        question: "Is there anything to take to a doctor?",
        sections: &[
            sec("clinvar", "ClinVar screen and hits", &["clinvar"]),
            sec(
                "clinical-filter",
                "Clinical filter and secondary findings",
                &["clinical"],
            ),
            sec("cpsr", "Cancer predisposition", &["cpsr"]),
            sec("slivar", "Variant prioritization", &["slivar"]),
            sec("exomiser", "Symptom-driven ranking", &["exomiser"]),
            sec(
                "structural-variants",
                "Structural variants",
                &["manta", "delly", "cnvpytor", "sv_consensus", "annotsv"],
            ),
            sec("mei", "Mobile element insertions", &["mei"]),
            sec("expansions", "Repeat expansions", &["expansions"]),
            sec("vntyper", "MUC1 repeat (ADTKD)", &["vntyper"]),
            sec("mito", "Mitochondrial variants", &["mito"]),
            sec("lpa", "Lp(a)", &["lpa"]),
            sec("chip", "Clonal haematopoiesis", &["chip"]),
            sec(
                "mosaic",
                "Mosaic loss of chromosome Y or X",
                &["mosaic_sex_loss"],
            ),
            sec("mosaic-mca", "Mosaic chromosomal alterations", &["mosaic"]),
        ],
    },
    Chapter {
        id: "carrier",
        title: "Carrier",
        question: "What could I pass on to a child?",
        sections: &[
            sec("carrier-list", "Carrier list", &["carrier"]),
            sec("couple", "Couple carrier screen", &["couple"]),
            sec(
                "paralogs",
                "Paralog genes and SMN1",
                &["paralogs", "parascopy"],
            ),
        ],
    },
    Chapter {
        id: "medicines",
        title: "Medicines",
        question: "How might I respond to medicines?",
        sections: &[
            sec("pharmcat", "PharmCAT", &["pharmcat"]),
            sec("cpic", "CPIC recommendations", &["cpic"]),
            sec("cyp2d6", "pypgx, Cyrius and CYP2D6", &["pypgx", "cyrius"]),
            sec("hla", "HLA typing", &["hla"]),
            sec("locityper", "Locityper second opinion", &["locityper"]),
        ],
    },
    Chapter {
        id: "ancestry",
        title: "Ancestry",
        question: "Where do my ancestors come from?",
        sections: &[
            sec("ancestry-fine", "Fine-grained ancestry", &["ancestry_fine"]),
            sec("mt-haplogroup", "Mitochondrial haplogroup", &["haplogroup"]),
            sec(
                "y-chromosome",
                "Y chromosome",
                &["y_haplogroup", "y_tree", "y_str"],
            ),
            sec("archaic", "Archaic ancestry", &["archaic"]),
            sec("ancient", "Ancient populations", &["ancient"]),
            sec("roh", "Runs of homozygosity", &["roh"]),
        ],
    },
    Chapter {
        id: "traits",
        title: "Traits and scores",
        question: "What do studies say about my traits and risk scores?",
        sections: &[
            sec("traits-panel", "Traits", &["traits"]),
            sec("prs", "Polygenic risk scores", &["prs"]),
            sec("bloodtype", "Blood groups", &["bloodtype"]),
            sec("telomere", "Telomere content", &["telomere"]),
        ],
    },
    Chapter {
        id: "other",
        title: "Other and exports",
        question: "What else did the run find or write?",
        sections: &[
            sec("viral", "Viral DNA", &["viral"]),
            sec("consumer-export", "Raw-data export", &["consumer_export"]),
        ],
    },
];

/// A section's text block: its lines, empty when the section is left out.
pub type TextFn = fn(&Summary) -> Result<Vec<String>, Malformed>;
/// A section's HTML part: its cards, none when the section is left out.
pub type HtmlFn = fn(&Summary) -> Result<HtmlPart, Malformed>;

/// The writers of section `sid`, or `None` for a section not ported yet,
/// which both reports leave out. One arm per section, a blank line apart, so
/// the units that port the others each change their own line.
pub fn writer(sid: &str) -> Option<(TextFn, HtmlFn)> {
    match sid {
        "qc" => Some((sections::qc::text, sections::qc::html)),

        "callability" => None,

        "variants" => Some((sections::variants::text, sections::variants::html)),

        "clinvar" => None,

        "clinical-filter" => None,

        "cpsr" => None,

        "slivar" => None,

        "exomiser" => None,

        "structural-variants" => None,

        "mei" => None,

        "expansions" => None,

        "vntyper" => None,

        "mito" => Some((sections::mito::text, sections::mito::html)),

        "lpa" => None,

        "chip" => None,

        "mosaic" => None,

        "mosaic-mca" => None,

        "carrier-list" => None,

        "couple" => None,

        "paralogs" => None,

        "pharmcat" => Some((sections::pgx::pharmcat_text, sections::pgx::pharmcat_html)),

        "cpic" => Some((sections::pgx::cpic_text, sections::pgx::cpic_html)),

        "cyp2d6" => Some((sections::pgx::cyp2d6_text, sections::pgx::cyp2d6_html)),

        "hla" => Some((sections::hla::text, sections::hla::html)),

        "locityper" => None,

        "ancestry-fine" => None,

        "mt-haplogroup" => None,

        "y-chromosome" => None,

        "archaic" => None,

        "ancient" => None,

        "roh" => None,

        "traits-panel" => None,

        "prs" => None,

        "bloodtype" => None,

        "telomere" => None,

        "viral" => None,

        "consumer-export" => None,

        _ => None,
    }
}

/// The chapter after the results: what the report leaves out.
pub const LIMITS: (&str, &str, &str) = (
    "limits",
    "What this report leaves out",
    "Which steps did not run, and what does no step look at?",
);

/// The closing appendix: every gene and rsID the results name.
pub const TERMS: (&str, &str, &str) = (
    "terms",
    "Genes and rsIDs in these results",
    "Where in this report does a gene or an rsID come up?",
);

/// Every summary key og summary collect writes, with its title and step,
/// in summary.json's order.
pub fn summary_keys() -> impl Iterator<Item = &'static Spec> {
    GROUP_A.iter().chain(GROUP_B.iter()).chain(GROUP_C.iter())
}

/// The step that writes summary key `key`, or `"?"` for a key no collector
/// writes.
pub fn step_of(key: &str) -> &'static str {
    summary_keys()
        .find(|s| s.key == key)
        .map_or("?", |s| s.step)
}

/// The section that shows summary key `key`.
pub fn section_of(key: &str) -> Option<&'static Sec> {
    CHAPTERS
        .iter()
        .flat_map(|c| c.sections.iter())
        .find(|s| s.keys.contains(&key))
}

/// The section with anchor `sid`.
pub fn section(sid: &str) -> Option<&'static Sec> {
    CHAPTERS
        .iter()
        .flat_map(|c| c.sections.iter())
        .find(|s| s.id == sid)
}

/// What is wrong with a chapter table against the summary keys: a key no
/// section shows, a key two sections show, a key no collector writes.
pub fn key_problems(chapters: &[Chapter], keys: &[&str]) -> Vec<String> {
    let seen: Vec<&str> = chapters
        .iter()
        .flat_map(|c| c.sections.iter().flat_map(|s| s.keys.iter().copied()))
        .collect();
    let mut out: Vec<String> = keys
        .iter()
        .filter(|k| !seen.contains(k))
        .map(|k| format!("missing {k}"))
        .collect();
    let mut twice: Vec<&str> = seen
        .iter()
        .copied()
        .filter(|k| seen.iter().filter(|x| *x == k).count() > 1)
        .collect();
    twice.sort_unstable();
    twice.dedup();
    out.extend(twice.iter().map(|k| format!("twice {k}")));
    out.extend(
        seen.iter()
            .filter(|k| !keys.contains(k))
            .map(|k| format!("unknown {k}")),
    );
    out
}

/// A section's state from the states of its keys. Never green: a section
/// that ran is not a finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameState {
    Ok,
    Stale,
    Failed,
    Part,
    Incomplete,
    Missing,
}

impl FrameState {
    /// The fold of the key states: failed or unreadable keys make it
    /// failed, or part failed beside a key that ran; a stale key makes it
    /// stale; a key that ran beside one that did not makes it part not run.
    pub fn of(states: &[State]) -> FrameState {
        let has = |s: State| states.contains(&s);
        let ran = has(State::Ok) || has(State::Stale);
        if has(State::Failed) || has(State::Unreadable) {
            return if ran {
                FrameState::Part
            } else {
                FrameState::Failed
            };
        }
        if has(State::Stale) {
            return FrameState::Stale;
        }
        match (has(State::Ok), has(State::Missing)) {
            (true, true) => FrameState::Incomplete,
            (true, false) => FrameState::Ok,
            (false, _) => FrameState::Missing,
        }
    }

    /// The class of the contents dot.
    pub fn class(self) -> &'static str {
        match self {
            FrameState::Ok => "ok",
            FrameState::Stale => "stale",
            FrameState::Failed => "failed",
            FrameState::Part => "part",
            FrameState::Incomplete => "incomplete",
            FrameState::Missing => "missing",
        }
    }

    /// The state badge under the first card title.
    pub fn badge(self) -> (&'static str, Colour) {
        match self {
            FrameState::Ok => ("Complete", Colour::State),
            FrameState::Stale => ("Stale", Colour::Yellow),
            FrameState::Failed => ("Failed", Colour::Red),
            FrameState::Part => ("Part failed", Colour::Red),
            FrameState::Incomplete => ("Part not run", Colour::Yellow),
            FrameState::Missing => ("Not run", Colour::Gray),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_state_folds_as_the_python_did() {
        use State::*;
        let got: Vec<FrameState> = [
            (Ok, Missing),
            (Stale, Missing),
            (Ok, Ok),
            (Missing, Missing),
            (Ok, Failed),
            (Unreadable, Missing),
            (Other, Ok),
        ]
        .iter()
        .map(|(a, b)| FrameState::of(&[*a, *b]))
        .collect();
        assert_eq!(
            got,
            [
                FrameState::Incomplete,
                FrameState::Stale,
                FrameState::Ok,
                FrameState::Missing,
                FrameState::Part,
                FrameState::Failed,
                FrameState::Ok,
            ]
        );
    }
}
