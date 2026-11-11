//! The processes of a module file (`modules/local/*/main.nf`): each one's
//! name, `publishDir` statements and `output:` paths, read from the
//! process's head (the part before its `script:`, `shell:`, `exec:` or
//! `stub:` section) with comments blanked first.

use crate::groovy::{closing_brace, strip_comments};
use regex::Regex;
use std::sync::LazyLock;

/// One `publishDir` statement, as written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublishDir {
    /// The first double-quoted string after `publishDir` (`""` when none).
    pub dir: String,
    /// The `pattern:` glob, when there is one.
    pub pattern: Option<String>,
    /// The body of the `saveAs:` closure, up to its own closing brace.
    pub save_as: Option<String>,
}

/// One process of a module file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Process {
    pub name: String,
    pub publish: Vec<PublishDir>,
    /// The quoted `path(...)` globs of its `output:` block, `versions.yml`
    /// left out.
    pub outputs: Vec<String>,
}

static PROCESS: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^process\s+(\w+)\s*\{"));
static SECTION: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^\s*(?:script|shell|exec|stub):"));
static PUBLISH: LazyLock<Regex> = LazyLock::new(|| re(r"^\s*publishDir\b"));
static DIR: LazyLock<Regex> = LazyLock::new(|| re(r#"publishDir\s+\{?\s*"([^"]+)""#));
static PATTERN: LazyLock<Regex> = LazyLock::new(|| re(r#"pattern:\s*["']([^"']+)["']"#));
static SAVE_AS: LazyLock<Regex> = LazyLock::new(|| re(r"saveAs:\s*\{"));
static OUTPUT: LazyLock<Regex> =
    LazyLock::new(|| re(r"(?ms)^\s*output:\s*$(.*?)(?:^\s*(?:when|topic):|\z)"));
static PATH: LazyLock<Regex> =
    LazyLock::new(|| re(r#"path\s*\(?\s*"([^"]+)"|path\s*\(?\s*'([^']+)'"#));

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

/// The body of a statement's `saveAs` closure, up to its own closing brace
/// (an option after it is not part of it); the rest of the statement when
/// the closure never closes; `None` when there is no `saveAs`.
pub fn closure_body(stmt: &str) -> Option<String> {
    let m = SAVE_AS.find(stmt)?;
    let open = m.end() - 1;
    Some(match closing_brace(stmt, open) {
        Some(close) => stmt[m.end()..close].to_owned(),
        None => stmt[m.end()..].to_owned(),
    })
}

/// The processes of the module source `src`, in file order.
pub fn scan(src: &str) -> Vec<Process> {
    let text = strip_comments(src);
    let starts: Vec<(usize, usize, String)> = PROCESS
        .captures_iter(&text)
        .map(|c| {
            let m = c.get(0).expect("match");
            (m.start(), m.end() - 1, c[1].to_owned())
        })
        .collect();
    let mut out = Vec::new();
    for (k, (start, open, name)) in starts.iter().enumerate() {
        let next = starts.get(k + 1).map_or(text.len(), |s| s.0);
        let end = closing_brace(&text, *open).map_or(next, |c| c.min(next));
        let body = &text[*start..end];
        let head = SECTION.find(body).map_or(body, |m| &body[..m.start()]);
        out.push(Process {
            name: name.clone(),
            publish: publish_dirs(head),
            outputs: outputs(head),
        });
    }
    out
}

/// The `publishDir` statements of a process head: a statement goes on to
/// the next line while it ends in `,` or has more `{` than `}`.
fn publish_dirs(head: &str) -> Vec<PublishDir> {
    let lines: Vec<&str> = head.split('\n').collect();
    let mut out = Vec::new();
    for (j, line) in lines.iter().enumerate() {
        if !PUBLISH.is_match(line) {
            continue;
        }
        let mut stmt = (*line).to_owned();
        let mut k = j;
        while stmt.trim_end().ends_with(',')
            || stmt.matches('{').count() > stmt.matches('}').count()
        {
            k += 1;
            let Some(l) = lines.get(k) else { break };
            stmt.push('\n');
            stmt.push_str(l);
        }
        out.push(PublishDir {
            dir: DIR
                .captures(&stmt)
                .map(|c| c[1].to_owned())
                .unwrap_or_default(),
            pattern: PATTERN.captures(&stmt).map(|c| c[1].to_owned()),
            save_as: closure_body(&stmt),
        });
    }
    out
}

fn outputs(head: &str) -> Vec<String> {
    let Some(block) = OUTPUT.captures(head) else {
        return Vec::new();
    };
    PATH.captures_iter(&block[1])
        .map(|c| {
            c.get(1)
                .or(c.get(2))
                .expect("one branch")
                .as_str()
                .to_owned()
        })
        .filter(|g| g != "versions.yml")
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SRC: &str = r#"
process ONE {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/x", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? null : fn }, enabled: { true }
    publishDir { "${params.outdir}/y" }, pattern: "*.{tsv,txt}"
    // publishDir "commented/out"

    output:
    tuple val(meta), path("*_a.tsv"), emit: a
    path('b.txt'), optional: true
    // path("commented.tsv")
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    publishDir "not/a/statement"
    echo '}' > c_a.tsv
    """
}

process TWO {
    output:
    path "two.vcf.gz"
    script:
    "true"
}
"#;

    #[test]
    fn reads_publish_dirs_and_outputs() {
        let procs = scan(SRC);
        assert_eq!(procs.len(), 2);
        let one = &procs[0];
        assert_eq!(one.name, "ONE");
        assert_eq!(one.outputs, ["*_a.tsv", "b.txt"]);
        assert_eq!(one.publish.len(), 2);
        assert_eq!(one.publish[0].dir, "${params.outdir}/${meta.id}/x");
        assert_eq!(
            one.publish[0].save_as.as_deref(),
            Some(" fn -> fn.endsWith('.log') ? null : fn ")
        );
        assert_eq!(one.publish[1].dir, "${params.outdir}/y");
        assert_eq!(one.publish[1].pattern.as_deref(), Some("*.{tsv,txt}"));
        assert_eq!(one.publish[1].save_as, None);
        assert_eq!(procs[1].name, "TWO");
        assert_eq!(procs[1].outputs, ["two.vcf.gz"]);
        assert!(procs[1].publish.is_empty());
    }

    #[test]
    fn a_closure_ends_at_its_own_brace() {
        let stmt =
            "publishDir \"x\", saveAs: { fn -> fn.endsWith('}') ? null : fn }, enabled: { true }";
        let body = closure_body(stmt).unwrap();
        assert!(!body.contains("enabled"), "{body}");
        assert!(body.contains("'}'"));
        assert_eq!(closure_body("publishDir \"x\""), None);
    }
}
