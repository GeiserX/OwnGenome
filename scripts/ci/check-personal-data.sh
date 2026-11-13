#!/usr/bin/env bash
# check-personal-data.sh: fail when a file holds something that looks like
# personal data: a home path, a private address, or a genotype.
#
# This repository is public and the pipeline runs on someone's own genome, so
# an output pasted into a doc, a test or a commit is a leak. Every tracked
# text file is scanned for:
#   home        a home path (/home/x, /Users/x, /mnt/user, C:\Users\x)
#   address     a private-range IPv4 address (10/8, 172.16/12, 192.168/16,
#               100.64/10) or a tailnet name (*.ts.net)
#   diplotype   a star-allele diplotype (*1/*4) on a line with a sample label
#               (sample, patient, subject, proband, NA12878, HG002, ...)
#   rsid        an rsID followed by a genotype (rs4244285 AG, rs123: A/G,
#               rs123 (C;T), or a chip export row rs123<TAB>1<TAB>123<TAB>AG)
#   repeat      a repeat-count genotype (17/19 on a line with REPCN, or after
#               an STR gene such as HTT or FMR1)
# Documented public examples (fixture values, the patterns themselves in the
# contributor guide) are listed in ALLOW below, by file and text; public or
# invented test data files are listed whole in SKIP. An entry
# covers only its own text: the rest of that line, and every other line of
# the file, is still checked. This script holds the list, so every listed
# text is allowed here as well.
#
# Usage:
#   scripts/ci/check-personal-data.sh             every tracked text file
#   scripts/ci/check-personal-data.sh FILE...     only these files
#   scripts/ci/check-personal-data.sh --cached FILE...
#                                                 the staged content of these
#                                                 files (git index), not the
#                                                 working tree
#   scripts/ci/check-personal-data.sh --self-test plant one line of each kind
#                                                 and require each to be caught
#
# Your own identifiers (names, sample ids, hostnames) do not belong in this
# public file. Put one Python-subset regular expression per text-file line
# outside the repository and point PERSONAL_DATA_PATTERNS at it. Each match
# is reported as kind "own". CRLF, CR and LF separate pattern-file lines.
#
# Supported: literals, dot, character classes and ranges (including negated
# classes), \b \B \w \W \s \S \d \D, ^ $ \A \Z, alternation, capturing,
# named (?P<name>...) and noncapturing groups, greedy and lazy ? * + {m,n},
# and global or scoped i/m/s/x flags. Hex \xHH, \uHHHH and \UHHHHHHHH and
# standard control escapes are supported. Unicode word characters, spaces
# and simple case folding follow Python, including dotted and dotless I.
# Every quantified operand must be nonnullable: it must consume a character
# on every successful match. This includes groups and alternatives under
# greedy or lazy ?, *, + and braces. For example, (a?)+ and (a*|b){1,3}
# are refused. Nullability propagates through flags, groups and alternatives.
# Empty matches outside such operands follow Python 3.14+ finditer. A
# nonempty match can follow an empty match at the same position. Captures never narrow the span checked against ALLOW.
#
# Backreferences, look-around in custom patterns, conditionals, atomic or
# possessive groups, ASCII/locale flags, named Unicode and octal escapes,
# nested/POSIX classes and other engine extensions are unsupported. Patterns
# over 16 KiB or with more than 64 nested groups are rejected. Unsupported
# syntax, invalid patterns and a matcher resource-limit error stop the scan
# with exit 2 and a diagnostic that omits the pattern and scanned values.
# A path with invalid UTF-8 also stops the scan before SKIP, path conversion
# or staged-blob requests. Symlink targets are also checked before SKIP:
# every working-copy hop (at most 40), staged targets from index mode 120000.
# Version-dependent filenames fail before ALLOW/SKIP with a generic error.
# Cached entries must be regular blobs or symlinks. Gitlinks and unreadable
# indexed objects fail with exit 2 naming only the path, unless SKIP excludes it.
# File contents still use replacement decoding.
#
# Unicode guard: UCD 13, 15, 16 and 17 disagree on newer script letters.
# The scan stops with redacted exit 2 at file:line for version-dependent
# text, filenames, custom patterns (including escaped literals), or ALLOW/SKIP rows.
# ALLOW cannot exempt such text. Use an ASCII SKIP glob to exempt its file.
# A letter from a script newer than Unicode 13 therefore stops the scan
# until its file is skipped. Every astral scalar with a case mapping or used
# as a mapping target in UCD 13/15/16/17 is also refused, regardless of the
# pattern shape. Text in Deseret, Osage, Old Hungarian, Warang Citi,
# Medefaidrin, Adlam and any other astral cased script stops the scan until
# its file is SKIPped. Newer plain symbols, emoji and uncased CJK still scan.
# A custom class range containing a property-guarded scalar is refused, even
# when its endpoints are stable. ASCII and stable Latin-1 ranges still work.
#
# The generated conservative set includes every scalar unassigned in UCD
# 17; post-13 letters, marks, numbers, connectors, spaces and case mappings;
# and older scalars whose word/digit/space membership or case maps change
# in UCD 13/15/16/17. Full case folds avoid rejecting six stable old scalars.
# BMP regex case differences apply only to patterns. Under active global or
# scoped i flags, a literal or a class range naming any member of a differing
# fold class is refused. This includes escaped members, positive and negated
# classes, and range interiors. Disabling i in that scope permits stable
# case-sensitive names and ranges. Ordinary Russian, Ukrainian and Serbian
# text, filenames and ALLOW/SKIP rows pass unless the property guard covers
# them. Dot and category escapes do not name a member of a fold class.
# Real fullmatch calls discover differing pairs. For each BMP pair, the
# converse verifier proves that property-guarded endpoints stop the scan.
# For admitted endpoints, it compares finditer spans using unrelated literals,
# categories, negated classes, ranges, alternations, scopes and quantifiers.
# It also discovers third case partners on every interpreter and checks their
# shapes. Any differing accepted shape on admitted text fails verification.
# Astral mappings come from UnicodeData, SpecialCasing
# and CaseFolding, all SHA-256-pinned, with no interpreter test needed.
# The 248 newer scalars with stable effective behavior remain conservatively
# blocked. The six stable scalars assigned by Unicode 13 pass.
# Text validation uses one binary search per non-ASCII character. Pattern
# compilation also checks a small fixed BMP fold-class set under i flags.
# That set is initialized once from the embedded pairs and engine fold table.
# This guard uses no per-file index or process. Cached symlink target names,
# including directory components, receive the same check before SKIP.
# Regenerate from SHA-256-pinned official data with the Rust generator:
# tests/parity/gen/personal-unicode-generate.sh SCRATCH_DIR
# Guard regenerates byte for byte and compares all available Python oracles,
# explicitly installing 3.9, 3.12 and 3.15. Match iteration follows 3.14+.
# Regenerate the regex table via the existing Rust differential test:
# PERSONAL_DATA_UPDATE_REGEX_CASES=1 cargo test -p og-dev round5_regex_case_table_matches_interpreters
# Supply all supported Python interpreters on PATH or PERSONAL_DATA_PYTHONS
# (colon-separated paths). Ordinary tests verify every available interpreter
# and reject any differing pair missing from the checked-in table. Every
# explicitly supplied path must run successfully. CI requires
# 3.9, 3.12 and 3.15 through PERSONAL_DATA_REQUIRE_PYTHONS. Missing explicit
# paths and failing interpreters are negative controls for both verifiers.
#
# The scan is og dev scan personal-data (rust/crates/og-dev); bin/og runs
# the build scripts/setup.sh --og puts in bin/.og for this machine, and OG
# names another og binary.
#
# As a pre-commit hook, checking the staged files:
#   cat > .git/hooks/pre-commit <<'HOOK'
#   #!/bin/sh
#   export PERSONAL_DATA_PATTERNS="$HOME/.config/pgp/personal-patterns"   # optional
#   git diff --cached --name-only --diff-filter=ACMR -z |
#     xargs -0 sh -c '[ "$#" -eq 0 ] || exec scripts/ci/check-personal-data.sh --cached "$@"' _
#   HOOK
#   chmod +x .git/hooks/pre-commit
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
# og runs the scan (bin/og picks the build for this machine); OG overrides it.
OG=${OG:-${ROOT}/bin/og}

# file (exact path or glob)   text on the line that makes it a documented example
# shellcheck disable=SC2016  # the entries are text, nothing expands
ALLOW='
CLAUDE.md                                     /mnt/user\|/Users/\|/home/
scripts/ci/check-personal-data.sh             (/home/x, /Users/x, /mnt/user, C:\Users\x)
scripts/ci/check-personal-data.sh             (*1/*4) on a line with a sample label
scripts/ci/check-personal-data.sh             (rs4244285 AG, rs123: A/G,
scripts/ci/check-personal-data.sh             rs123 (C;T), or a chip export row
scripts/ci/check-personal-data.sh             (17/19 on a line with REPCN
CONTRIBUTING.md                               /mnt/user\|internal-host\|/home/
CONTRIBUTING.md                               (`/mnt/user/`, `/home/username/`, etc.)
modules/local/cyrius/main.nf                  ${prefix}\\t*1/*1\\tPASS
tests/fake-docker/reports-stale-step.sh       sample1\t*1/*1\tPASS
tests/test_collect_summary.py                 S\t*1/*4\tPASS
tests/test_collect_summary.py                 strec("HTT", "17/30"
tests/test_collect_summary.py                 strec("C9ORF72", "2/8"
tests/test_collect_summary.py                 strec("ATXN1", "29/31"
tests/test_collect_summary.py                 strec("DMPK", "5/40"
tests/test_collect_summary.py                 strec("RFC1", "100/150"
tests/smoke/commands.tsv                      HTT at 19/45 repeats
tests/smoke/commands.tsv                      ATXN3 at 22/24 repeats
'

# Whole files that are not scanned (glob), with the reason. Keep this list
# short: an entry hides every line of the file.
SKIP='
tests/fixtures/*                public or invented test data; each fixture directory says where it came from
tests/smoke/*.vcf.in            invented tool inputs for the image smoke tests
tests/parity/inputs/pgx-refblock-fill/*  built from the public PharmCAT data of tests/fixtures/pharmcat
tests/parity/golden/pgx-refblock-fill/*  the output og writes on those inputs
'

# scan [FILE...]: og dev scan personal-data (rust/crates/og-dev) with the
# lists above; with no FILE, every file git ls-files lists under ROOT. A
# scan that read no text file fails: it checked nothing.
scan() {
  local cached=()
  [ -z "${CACHED:-}" ] || cached=(--cached)
  local patterns=()
  [ -z "${PERSONAL_DATA_PATTERNS:-}" ] || patterns=(--patterns "$PERSONAL_DATA_PATTERNS")
  "$OG" dev scan personal-data --root "$ROOT" --allow "$ALLOW" --skip "$SKIP" \
    ${patterns[@]+"${patterns[@]}"} ${cached[@]+"${cached[@]}"} -- "$@"
}

self_test() {
  local tmp fail=0 rc out k
  tmp=$(mktemp -d)
  # shellcheck disable=SC2064
  trap "rm -rf '$tmp'" EXIT
  # One planted line per kind, each in its own file. All values are made up.
  # plant drops every "~", which keeps this file itself free of the patterns.
  plant() { printf '%s\n' "${2//\~/}" > "${tmp}/$1"; }
  plant home.md 'cp results.html /home/~alice/genome/'
  plant address.md 'scp out.vcf.gz root@192.~168.7.20:/data/'
  plant tailnet.md 'ssh nas.tail1234.ts.~net'
  plant diplotype.md '| sample7 | CYP2D6 | *1/~*4 |'
  plant label.md 'sample_01 CYP2D6 *1/~*4'
  plant rsid.md 'rs4244285~ AG'
  printf 'rs4477212\t1\t82154\tAA\n' > "${tmp}/chip.txt"
  plant repeat.md 'HTT CAG repeats: 17/~19'
  plant repcn.vcf 'chr4 3074877 . C <STR17> . PASS . GT:REPCN 0/~1:17/~19'
  for k in home address tailnet diplotype label rsid chip repeat repcn; do
    f=$(ls "${tmp}/${k}".*)
    rc=0; out=$(scan "$f" 2>&1) || rc=$?
    if [ "$rc" -ne 1 ] || ! grep -q '^FAIL' <<<"$out"; then
      echo "self-test: planted ${k} line was not caught (exit ${rc}):"; printf '%s\n' "$out"; fail=1
    else
      echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
    fi
  done
  # Lines that must pass: doc text without a sample, a version, a public URL.
  # shellcheck disable=SC2016  # the backticks are text
  printf '%s\n' '- `*1/*4` -- Intermediate metabolizer' 'samtools 1.21, GATK 4.6.2.0' \
    'https://example.org/home/page' 'CYP2C19*2 is rs4244285 (G>A)' > "${tmp}/clean.md"
  rc=0; out=$(scan "${tmp}/clean.md" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: clean lines were reported:"; printf '%s\n' "$out"; fail=1; }
  # An allowed example covers its own text only: alone it passes, with a home
  # path appended to the same line it is reported.
  mkdir -p "${tmp}/tree"
  # shellcheck disable=SC2016  # the backticks are text
  plant tree/CONTRIBUTING.md 'Personal file paths (`/mnt/user/`, `/home/username/`, etc.)'
  rc=0; out=$(ROOT="${tmp}/tree" scan "${tmp}/tree/CONTRIBUTING.md" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: an allowed example was reported:"; printf '%s\n' "$out"; fail=1; }
  # shellcheck disable=SC2016  # the backticks are text
  plant tree/CONTRIBUTING.md 'Personal file paths (`/mnt/user/`, `/home/username/`, etc.) see /home/~bob/x'
  rc=0; out=$(ROOT="${tmp}/tree" scan "${tmp}/tree/CONTRIBUTING.md" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL home' <<<"$out"; then
    echo "self-test: a home path after an allowed example was not caught:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
  # --cached reads the index: a staged leak is caught after the working copy
  # drops it.
  git -C "${tmp}/tree" init -q
  plant tree/staged.md 'cp results.html /home/~alice/genome/'
  git -C "${tmp}/tree" add staged.md
  plant tree/staged.md 'cp results.html genome/'
  rc=0; out=$(cd "${tmp}/tree" && ROOT="${tmp}/tree" CACHED=1 scan staged.md 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL home' <<<"$out"; then
    echo "self-test: --cached did not read the staged content:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
  # A scan that read no text file checked nothing: it fails. Here the only
  # file is binary.
  printf 'a\0b\n' > "${tmp}/binary.bin"
  rc=0; out=$(scan "${tmp}/binary.bin" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL: no text file was scanned' <<<"$out"; then
    echo "self-test: a scan of no text file passed (exit ${rc}):"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
  # ROOT reached through a symbolic link still places the file below it, so
  # the file's allowed example applies.
  ln -s "${tmp}/tree" "${tmp}/link"
  # shellcheck disable=SC2016  # the backticks are text
  plant tree/CONTRIBUTING.md 'Personal file paths (`/mnt/user/`, `/home/username/`, etc.)'
  rc=0; out=$(ROOT="${tmp}/link" scan "${tmp}/tree/CONTRIBUTING.md" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: an allowed example was reported under a linked ROOT:"; printf '%s\n' "$out"; fail=1; }
  # A pattern of your own, from PERSONAL_DATA_PATTERNS.
  printf '%s\n' 'Report for Jane Example' > "${tmp}/own.md"
  printf '%s\n' '\bJane Example\b' > "${tmp}/patterns"
  rc=0; out=$(PERSONAL_DATA_PATTERNS="${tmp}/patterns" scan "${tmp}/own.md" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q '^FAIL own' <<<"$out"; then
    echo "self-test: a PERSONAL_DATA_PATTERNS match was not caught:"; printf '%s\n' "$out"; fail=1
  else
    echo "self-test: caught: $(grep '^FAIL' <<<"$out")"
  fi
  [ "$fail" -eq 0 ] && echo "self-test: OK"
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  --cached) shift; CACHED=1 scan "$@" ;;
  -*) echo "usage: $0 [FILE... | --self-test]" >&2; exit 2 ;;
  *) scan "$@" ;;
esac
