//! A shell lexer just good enough to find the commands a bash script runs:
//! the words of every simple command and where the command starts.
//!
//! It knows quotes, `$(…)` and `<(…)`/`>(…)` (their commands are lexed too,
//! and the substitution itself becomes the word `$(...)`), `${…}`,
//! comments, here-documents (their bodies are skipped, `<<-` strips the
//! leading tabs), redirections (the operator and its target are dropped)
//! and the separators `;`, `|`, `&`, `(`, `)`, `&&`, `||`, `;;` and `|&`.
//! It is not a shell: `$((…))` arithmetic splits a command at its `(`, and
//! a `}` inside quotes in `${…}` closes the expansion early. Neither occurs
//! in a `docker run` or `run_in` call of this repository.
//!
//! The lexer works on bytes; every offset is a byte offset into the text it
//! was given.

/// One simple command: its words, quotes removed, and the byte offset of
/// its first word.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Command {
    pub words: Vec<String>,
    pub offset: usize,
}

/// Words that start a compound command or prefix a simple one; the command
/// word comes after them.
const KEYWORDS: [&str; 17] = [
    "if", "then", "elif", "else", "fi", "do", "done", "while", "until", "{", "}", "!", "time",
    "function", "exec", "[[", "]]",
];

/// Words after which the rest of the command is not a command (`for x in
/// …`, `case … in`).
const STOP: [&str; 5] = ["for", "select", "case", "in", "esac"];

/// A keyword or prefix the command word comes after.
pub fn is_keyword(w: &str) -> bool {
    KEYWORDS.contains(&w)
}

/// `NAME=…`, `NAME+=…` or `NAME[i]=…`: an assignment before the command.
pub fn is_assignment(w: &str) -> bool {
    let b = w.as_bytes();
    if !b
        .first()
        .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'_')
    {
        return false;
    }
    let mut k = 1;
    while k < b.len() && (b[k].is_ascii_alphanumeric() || b[k] == b'_') {
        k += 1;
    }
    if b.get(k) == Some(&b'[') {
        match b[k..].iter().position(|&c| c == b']') {
            Some(close) => k += close + 1,
            None => return false,
        }
    }
    if b.get(k) == Some(&b'+') {
        k += 1;
    }
    b.get(k) == Some(&b'=')
}

/// The command words of `cmds`: for each command, the first word after the
/// keywords and assignments in front. A command whose first such word is a
/// path (it holds a `/`) runs that file, not a binary from `PATH`, and
/// yields nothing; so does one that starts with `for`, `case` and the like.
pub fn heads(cmds: &[Command]) -> Vec<String> {
    let mut out = Vec::new();
    for cmd in cmds {
        for w in &cmd.words {
            if STOP.contains(&w.as_str()) {
                break;
            }
            if is_keyword(w) || is_assignment(w) {
                continue;
            }
            if !w.contains('/') {
                out.push(w.clone());
            }
            break;
        }
    }
    out
}

/// The simple commands of `text`, in the order their ends are reached: a
/// command inside a substitution comes before the command around it.
pub fn lex(text: &str) -> Vec<Command> {
    let mut cmds = Vec::new();
    lex_into(text.as_bytes(), 0, &mut cmds);
    cmds
}

fn find(t: &[u8], c: u8, from: usize) -> Option<usize> {
    t.get(from..)?
        .iter()
        .position(|&x| x == c)
        .map(|k| from + k)
}

/// `t[a..b]`, empty when `b` is not past `a`.
fn slice(t: &[u8], a: usize, b: usize) -> &[u8] {
    if b > a && a <= t.len() {
        &t[a..b.min(t.len())]
    } else {
        &[]
    }
}

/// `t[j]` is a double quote: the index after its closing quote.
fn skip_dq(t: &[u8], mut j: usize) -> usize {
    j += 1;
    while j < t.len() && t[j] != b'"' {
        if t[j] == b'\\' {
            j += 2;
        } else if t[j..].starts_with(b"$(") {
            j = match_paren(t, j + 1) + 1;
        } else {
            j += 1;
        }
    }
    j + 1
}

/// `t[i]` is `(`: the index of its matching `)`, or the last index when it
/// has none.
pub(crate) fn match_paren(t: &[u8], i: usize) -> usize {
    let mut depth = 0i32;
    let mut j = i;
    while j < t.len() {
        match t[j] {
            b'\\' => {
                j += 2;
                continue;
            }
            b'\'' => {
                j = find(t, b'\'', j + 1).map_or(t.len(), |k| k + 1);
                continue;
            }
            b'"' => {
                j = skip_dq(t, j);
                continue;
            }
            b'(' => depth += 1,
            b')' => {
                depth -= 1;
                if depth == 0 {
                    return j;
                }
            }
            _ => {}
        }
        j += 1;
    }
    t.len().saturating_sub(1)
}

/// `t[i]` is `{`: the index of its matching `}`, or the last index when it
/// has none. Quotes are not looked at.
pub(crate) fn match_brace(t: &[u8], i: usize) -> usize {
    let mut depth = 0i32;
    for (j, &c) in t.iter().enumerate().skip(i) {
        if c == b'{' {
            depth += 1;
        } else if c == b'}' {
            depth -= 1;
            if depth == 0 {
                return j;
            }
        }
    }
    t.len().saturating_sub(1)
}

#[derive(Default)]
struct State {
    /// The words of the command being read, with their offsets.
    cur: Vec<(Vec<u8>, usize)>,
    /// The word being read, and its offset.
    word: Option<(Vec<u8>, usize)>,
    /// The next word to end is a redirection target and is dropped.
    redirect: bool,
    /// Here-documents whose bodies start at the next newline: delimiter
    /// and whether `<<-` strips leading tabs.
    heredocs: Vec<(Vec<u8>, bool)>,
}

impl State {
    fn add(&mut self, s: &[u8], at: usize) {
        self.word
            .get_or_insert_with(|| (Vec::new(), at))
            .0
            .extend_from_slice(s);
    }

    fn end_word(&mut self) {
        if let Some(w) = self.word.take() {
            if self.redirect {
                self.redirect = false;
            } else {
                self.cur.push(w);
            }
        }
    }

    fn end_cmd(&mut self, base: usize, cmds: &mut Vec<Command>) {
        self.end_word();
        if let Some(&(_, first)) = self.cur.first() {
            cmds.push(Command {
                words: self
                    .cur
                    .drain(..)
                    .map(|(w, _)| String::from_utf8_lossy(&w).into_owned())
                    .collect(),
                offset: base + first,
            });
        }
    }
}

/// A here-document line ends the body when it is the delimiter itself, as
/// bash reads it: the whole line, after the leading tabs for `<<-`.
fn ends_heredoc(line: &[u8], delim: &[u8], strip_tabs: bool) -> bool {
    let line = if strip_tabs {
        let k = line.iter().position(|&c| c != b'\t').unwrap_or(line.len());
        &line[k..]
    } else {
        line
    };
    line == delim
}

/// The delimiter word after `<<` or `<<-` (and blanks) at `t[j..]`: `NAME`,
/// `'NAME'` or `"NAME"`, with the index after it.
fn heredoc_word(t: &[u8], j: usize) -> Option<(Vec<u8>, usize)> {
    let quote = t.get(j).copied().filter(|&c| c == b'\'' || c == b'"');
    let start = j + usize::from(quote.is_some());
    let mut end = start;
    while end < t.len() && (t[end].is_ascii_alphanumeric() || t[end] == b'_') {
        end += 1;
    }
    if end == start {
        return None;
    }
    match quote {
        None => Some((t[start..end].to_vec(), end)),
        Some(q) if t.get(end) == Some(&q) => Some((t[start..end].to_vec(), end + 1)),
        Some(_) => None,
    }
}

fn lex_into(t: &[u8], base: usize, cmds: &mut Vec<Command>) {
    let n = t.len();
    let mut st = State::default();
    let mut i = 0;
    while i < n {
        let c = t[i];
        let rest = &t[i..];
        if c == b'\\' {
            if t.get(i + 1) == Some(&b'\n') {
                i += 2;
                continue;
            }
            st.add(slice(t, i + 1, i + 2), i);
            i += 2;
        } else if c == b'\'' {
            let k = find(t, b'\'', i + 1).unwrap_or(n);
            st.add(slice(t, i + 1, k), i);
            i = k + 1;
        } else if c == b'"' {
            let mut j = i + 1;
            let mut buf = Vec::new();
            while j < n && t[j] != b'"' {
                if t[j] == b'\\' && j + 1 < n && b"\"$\\`\n".contains(&t[j + 1]) {
                    if t[j + 1] != b'\n' {
                        buf.push(t[j + 1]);
                    }
                    j += 2;
                    continue;
                }
                if t[j..].starts_with(b"$(") && !t[j..].starts_with(b"$((") {
                    let k = match_paren(t, j + 1);
                    lex_into(slice(t, j + 2, k), base + j + 2, cmds);
                    buf.extend_from_slice(b"$(...)");
                    j = k + 1;
                    continue;
                }
                buf.push(t[j]);
                j += 1;
            }
            st.add(&buf, i);
            i = j + 1;
        } else if (rest.starts_with(b"$(") && !rest.starts_with(b"$(("))
            || ((c == b'<' || c == b'>') && t.get(i + 1) == Some(&b'('))
        {
            let k = match_paren(t, i + 1);
            lex_into(slice(t, i + 2, k), base + i + 2, cmds);
            st.add(b"$(...)", i);
            i = k + 1;
        } else if rest.starts_with(b"${") {
            let k = match_brace(t, i + 1);
            st.add(slice(t, i, k + 1), i);
            i = k + 1;
        } else if c == b'#' && st.word.is_none() {
            i = find(t, b'\n', i).unwrap_or(n);
        } else if c == b' ' || c == b'\t' {
            st.end_word();
            i += 1;
        } else if c == b'\n' {
            st.end_cmd(base, cmds);
            i += 1;
            for (delim, strip) in std::mem::take(&mut st.heredocs) {
                while i < n {
                    let k = find(t, b'\n', i).unwrap_or(n);
                    let line = &t[i..k];
                    i = k + 1;
                    if ends_heredoc(line, &delim, strip) {
                        break;
                    }
                }
            }
        } else if b";|&()".contains(&c) {
            if rest.starts_with(b"&>") {
                st.end_word();
                st.redirect = true;
                i += 2;
            } else if [&b"||"[..], b"&&", b";;", b"|&"]
                .iter()
                .any(|p| rest.starts_with(p))
            {
                st.end_cmd(base, cmds);
                i += 2;
            } else {
                st.end_cmd(base, cmds);
                i += 1;
            }
        } else if c == b'<' || c == b'>' {
            if rest.starts_with(b"<<<") {
                st.end_word();
                i += 3;
            } else if rest.starts_with(b"<<") {
                st.end_word();
                let mut j = i + 2;
                let mut strip = false;
                if t.get(j) == Some(&b'-') {
                    strip = true;
                    j += 1;
                }
                while j < n && (t[j] == b' ' || t[j] == b'\t') {
                    j += 1;
                }
                if let Some((delim, after)) = heredoc_word(t, j) {
                    st.heredocs.push((delim, strip));
                    j = after;
                }
                i = j;
            } else {
                // A file descriptor number before the operator is part of it.
                if st
                    .word
                    .as_ref()
                    .is_some_and(|(w, _)| !w.is_empty() && w.iter().all(u8::is_ascii_digit))
                {
                    st.word = None;
                }
                st.end_word();
                let mut j = i;
                while j < n && b"<>&|".contains(&t[j]) {
                    j += 1;
                }
                st.redirect = true;
                i = j;
            }
        } else {
            st.add(&[c], i);
            i += 1;
        }
    }
    st.end_cmd(base, cmds);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(text: &str) -> Vec<Vec<String>> {
        lex(text).into_iter().map(|c| c.words).collect()
    }

    fn heads_of(text: &str) -> Vec<String> {
        heads(&lex(text))
    }

    #[test]
    fn splits_commands_and_drops_redirections() {
        assert_eq!(
            words("a x | b 2>/dev/null && c >out; d &> log || e <in\n"),
            [vec!["a", "x"], vec!["b"], vec!["c"], vec!["d"], vec!["e"]]
        );
    }

    #[test]
    fn quotes_join_and_escape() {
        assert_eq!(
            words(r#"echo 'a b' "c \"d\" \$e" f\ g"#),
            [vec!["echo", "a b", "c \"d\" $e", "f g"]]
        );
    }

    #[test]
    fn substitutions_are_lexed_first() {
        let cmds = lex("OUT=$(time docker run img tabix z)\n");
        assert_eq!(
            cmds[0].words,
            ["time", "docker", "run", "img", "tabix", "z"]
        );
        assert_eq!(cmds[0].offset, 6);
        assert_eq!(cmds[1].words, ["OUT=$(...)"]);
        assert_eq!(
            heads_of("x=\"$(bgzip -c a)\"; diff <(tabix -l b) c"),
            ["bgzip", "tabix", "diff"]
        );
    }

    #[test]
    fn heads_skip_keywords_assignments_and_paths() {
        assert_eq!(heads_of("if ! bgzip x; then"), ["bgzip"]);
        assert_eq!(heads_of("time tabix y"), ["tabix"]);
        assert_eq!(heads_of("VAR=x A[1]=2 B+=3 samtools view"), ["samtools"]);
        assert_eq!(heads_of("\"$d/bgzip\" -c x"), Vec::<String>::new());
        assert_eq!(heads_of("for f in a b; do tabix $f; done"), ["tabix"]);
    }

    #[test]
    fn comments_and_brace_expansions() {
        assert_eq!(words("a ${x:-b c} # bgzip\n"), [vec!["a", "${x:-b c}"]]);
        assert_eq!(words("a#b"), [vec!["a#b"]]);
    }

    #[test]
    fn here_documents_are_skipped() {
        let text = "cat <<EOF > f\nbgzip inside\nEOF\ntabix after\n";
        assert_eq!(heads_of(text), ["cat", "tabix"]);
        let tabs = "cat <<-'END'\n\tbgzip inside\n\tEND\ntabix after\n";
        assert_eq!(heads_of(tabs), ["cat", "tabix"]);
    }

    #[test]
    fn a_here_document_ends_only_on_its_exact_line() {
        // Bash does not end a <<EOF body on "  EOF" (leading blanks) or on
        // "EOF " (a trailing blank): the body goes on to the real "EOF".
        let text = "cat <<EOF\n  EOF\nbgzip inside\nEOF \nEOF\ntabix after\n";
        assert_eq!(heads_of(text), ["cat", "tabix"]);
    }

    #[test]
    fn offsets_are_byte_offsets() {
        let text = "é=1\n  bgzip x\n";
        let cmds = lex(text);
        assert_eq!(cmds[1].words, ["bgzip", "x"]);
        assert_eq!(&text[cmds[1].offset..cmds[1].offset + 5], "bgzip");
    }

    #[test]
    fn unclosed_constructs_do_not_panic() {
        for t in [
            "$(", "\"$(", "${", "'", "\"", "<(", "cat <<", "a \\", "$(a 'b",
        ] {
            let _ = lex(t);
        }
    }

    #[test]
    fn assignment_shapes() {
        assert!(is_assignment("A=1"));
        assert!(is_assignment("_a[x y]+=1"));
        assert!(!is_assignment("1A=1"));
        assert!(!is_assignment("--opt=1"));
        assert!(!is_assignment("A[=1"));
    }
}
