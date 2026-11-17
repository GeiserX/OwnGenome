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
#   1. a planted page with such a row fails, and the scanner names its line;
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
# a backslash and a pipe; returns 1 when there is one. A table row starts with
# a pipe. Lines inside ``` or ~~~ fences are skipped. A code span opens with a
# run of N backticks and closes at the next run of exactly N; a run with no
# partner is plain text.
scan() {
  awk '
    function bad(line,   rest, tmp, n, off, found, code) {
      rest = line
      while (match(rest, /`+/)) {
        n = RLENGTH
        rest = substr(rest, RSTART + n)
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
    FNR == 1 { fence = 0 }
    /^[ \t]*(```|~~~)/ { fence = !fence; next }
    fence { next }
    /^[ \t]*\|/ && bad($0) { print FILENAME ":" FNR ": " $0; hits++ }
    END { exit hits > 0 }
  ' "$@"
}

# 1. A planted row like the old Nextflow line fails, and the output names it.
cat > "${WORK}/bad.md" <<'MD'
# Planted

| Software | Install |
|---|---|
| Docker | from docs.docker.com |
| Nextflow | `curl -s https://get.nextflow.io \| NXF_VER=1.0.0 bash` |
MD
if out=$(scan "${WORK}/bad.md"); then
  fail "a planted table row with an escaped pipe in a code span passed"
elif grep -qF "bad.md:6: | Nextflow |" <<<"$out"; then
  pass "a planted table row with an escaped pipe in a code span fails and is named"
else
  fail "the planted row failed but its line was not named: ${out}"
fi

# Double backticks open a span too.
printf '| a | ``x \\| y`` |\n' > "${WORK}/double.md"
if scan "${WORK}/double.md" >/dev/null; then
  fail "an escaped pipe in a double-backtick span passed"
else
  pass "an escaped pipe in a double-backtick span fails"
fi

# 2. Pipes where they render right pass.
cat > "${WORK}/good.md" <<'MD'
# Planted

Prose may say `a \| b` outside a table.

| Column | Meaning |
|---|---|
| `readings` | `KEY=text`, separated by a vertical bar |
| choice | `left` \| `right` |
| stray | a lone `` run, then \| and `code` |

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

# 3. The real docs.
mapfile -t pages < <(find "$DOCS" -name '*.md' -type f | sort)
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
