#!/usr/bin/env bash
# A command inside a code span in a docs table must not hold an escaped pipe.
#
# GitHub needs `\|` to keep a pipe inside a table cell, but the docs site
# prints a code span as written, so the reader sees and copies the backslash:
# `curl -s https://get.nextflow.io \| NXF_VER=... bash` is not a pipeline.
# A command that needs a pipe goes in a code block under the table, and a
# table cell names a separator in words.
#
# Checks:
#   1. planted pages with such a row fail, and the scanner names the line:
#      a bordered table, a double-backtick span, a table without border
#      pipes, a span after an escaped backtick, and a table after a fenced
#      block that holds a shorter fence line;
#   2. a planted page with a pipe only where it renders right passes: a
#      fenced block, prose, a table cell outside a code span, a code span in
#      prose, and backticks that do not pair;
#   3. every docs/**/*.md page passes.
#
# Run: bash tests/test_docs_table_pipes.sh [DIR]   (DIR defaults to docs/)
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
DOCS=${1:-"${REPO}/docs"}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# scan FILE...: prints FILE:LINE: TEXT for each table row whose code span holds
# a backslash and a pipe; returns 1 when there is one. A table starts at a line
# with a pipe followed by a separator line (`|---|---|` or `--- | ---`, with or
# without border pipes) and ends at a blank line; its header and every row are
# checked. A fenced block opens with 3 or more backticks or tildes and closes
# only at a run of the same character at least as long with nothing after it;
# lines inside are skipped. A code span opens with a run of N backticks and
# closes at the next run of exactly N; a run with no partner is plain text, and
# a backtick after a backslash is a literal, not an opener.
scan() {
  awk '
    function runlen(s, c,   i) { i = 0; while (substr(s, i + 1, 1) == c) i++; return i }
    function bad(line,   rest, tmp, n, s, k, off, found, code) {
      rest = line
      while (match(rest, /`+/)) {
        n = RLENGTH; s = RSTART
        k = 0; while (s - k > 1 && substr(rest, s - k - 1, 1) == "\\") k++
        if (k % 2) { s++; n-- }
        rest = substr(rest, s + n)
        if (n == 0) continue
        tmp = rest; off = 0; found = 0
        while (match(tmp, /`+/)) {
          if (RLENGTH == n) {
            code = substr(rest, 1, off + RSTART - 1)
            rest = substr(tmp, RSTART + RLENGTH)
            found = 1
            break
          }
          off += RSTART + RLENGTH - 1
          tmp = substr(tmp, RSTART + RLENGTH)
        }
        if (found && index(code, "\\|") > 0) return 1
      }
      return 0
    }
    function hit(n, text) { print FILENAME ":" n ": " text; hits++ }
    FNR == 1 { fence = 0; table = 0; prev = ""; prevn = 0 }
    fence {
      t = $0; sub(/^[ \t]*/, "", t); n = runlen(t, fch)
      if (n >= flen && substr(t, n + 1) ~ /^[ \t]*$/) fence = 0
      next
    }
    /^[ \t]*(```|~~~)/ {
      t = $0; sub(/^[ \t]*/, "", t); fch = substr(t, 1, 1); flen = runlen(t, fch)
      fence = 1; table = 0; prev = ""
      next
    }
    /^[ \t]*$/ { table = 0; prev = ""; next }
    table { if (bad($0)) hit(FNR, $0); next }
    /\|/ && /^[ \t]*\|?[ \t]*:?-+:?[ \t]*(\|[ \t]*:?-+:?[ \t]*)*\|?[ \t]*$/ && index(prev, "|") > 0 {
      table = 1
      if (bad(prev)) hit(prevn, prev)
      next
    }
    { prev = $0; prevn = FNR }
    END { exit hits > 0 }
  ' "$@"
}

# expect_hit NAME LINE WHAT: the planted page NAME must fail, and the scanner
# must name exactly one row, at line LINE. A scanner error that names no row
# is a failure too.
expect_hit() {
  local name=$1 line=$2 what=$3 out rc=0
  out=$(scan "${WORK}/${name}" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ]; then
    fail "${what} passed"
  elif [ "$(grep -c "^${WORK}/${name}:" <<<"$out")" -eq 1 ] && grep -q "^${WORK}/${name}:${line}: " <<<"$out"; then
    pass "${what} fails and is named"
  else
    fail "${what} failed, but the scanner did not name only line ${line}: ${out}"
  fi
}

# 1. Planted rows that must fail.
# A row like the old Nextflow line.
cat > "${WORK}/bad.md" <<'MD'
# Planted

| Software | Install |
|---|---|
| Docker | from docs.docker.com |
| Nextflow | `curl -s https://get.nextflow.io \| NXF_VER=1.0.0 bash` |
MD
expect_hit bad.md 6 "a planted table row with an escaped pipe in a code span"

# Double backticks open a span too.
cat > "${WORK}/double.md" <<'MD'
| a | b |
|---|---|
| a | ``x \| y`` |
MD
expect_hit double.md 3 "an escaped pipe in a double-backtick span"

# A table without border pipes is a table too.
cat > "${WORK}/borderless.md" <<'MD'
Tool | Install
--- | ---
Docker | from docs.docker.com
Nextflow | `curl -s https://get.nextflow.io \| bash`
MD
expect_hit borderless.md 4 "an escaped pipe in a code span of a table without border pipes"

# An escaped backtick is a literal, so the span starts at the next backtick.
cat > "${WORK}/escaped.md" <<'MD'
| a | b |
|---|---|
| a | \`example `x \| y` |
MD
expect_hit escaped.md 3 "an escaped pipe in a code span after an escaped backtick"

# A shorter fence line inside a four-backtick block does not close it, so the
# row inside stays skipped and the table after the block is checked.
cat > "${WORK}/fence.md" <<'MD'
````markdown
```
| `inside the block \| here` |
````

| a | b |
|---|---|
| a | `x \| y` |
MD
expect_hit fence.md 8 "an escaped pipe in a table after a block with a shorter fence line"

# 2. Pipes where they render right pass.
cat > "${WORK}/good.md" <<'MD'
# Planted

Prose may say `a \| b` outside a table.

| Column | Meaning |
|---|---|
| `readings` | `KEY=text`, separated by a vertical bar |
| choice | `left` \| `right` |
| stray | a lone `` run, then \| and `code` |

A line that starts with a pipe but has no separator under it is not a table:
| `not a table \| row` |

```bash
dmesg | grep -i "oom\|killed"
| `not a table \| row` |
```

~~~
| `still fenced \| here` |
~~~
MD
if out=$(scan "${WORK}/good.md"); then
  pass "pipes in fences, prose, cells outside code spans and unpaired backticks pass"
else
  fail "a page with no bad row failed: ${out}"
fi

# 3. The real docs. A find that cannot read part of the tree fails the test,
# so a partial list never passes as a clean scan.
if ! find "$DOCS" -name '*.md' -type f > "${WORK}/pages.raw"; then
  fail "could not list the Markdown pages under ${DOCS}"
fi
pages=()
while IFS= read -r page; do pages+=("$page"); done < <(sort "${WORK}/pages.raw")
if [ "${#pages[@]}" -eq 0 ]; then
  fail "no Markdown pages under ${DOCS}"
elif out=$(scan "${pages[@]}"); then
  pass "no docs table row has an escaped pipe in a code span (${#pages[@]} pages)"
else
  echo "$out"
  fail "docs table rows with an escaped pipe in a code span; move the command into a code block under the table"
fi

if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All checks passed."
