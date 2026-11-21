#!/usr/bin/env bash
# A new user must see what to install before the first command to copy, and
# one ordered install path for their platform.
#
# README.md:
#   - a "Before you begin" list (Docker, bash, Java, Nextflow, RAM) comes
#     before the first fenced code block;
#   - its RAM line says a 16 GB machine needs swap for step 25 (25-prs #memory);
#   - the first code block clones the repository and checks out the newest
#     release tag, so setup downloads og instead of building it. The tag is
#     picked with exactly git tag --list 'v[0-9]*' --sort=-v:refname, so a
#     fixture-* or e2e-data-* tag never wins.
# docs/getting-started.md:
#   - the RAM row of the hardware table names swap for step 25 (25-prs.md#memory);
#   - "### Get the code" clones and checks out the release tag the same way;
#   - four sections, "### Linux", "### Intel Mac", "### Apple Silicon Mac" and
#     "### Windows with WSL2", each installing Java and the pinned Nextflow,
#     with Docker named; each section's last code block is the verification
#     block (docker info, bash --version, java -version, nextflow -version)
#     and a "Next" line follows it;
#   - inside each section the steps run in an order that works on a new
#     machine: the Java install line comes before the first NXF_VER= line
#     (the Nextflow installer needs Java); in the two Mac sections the
#     Homebrew install line comes before the first brew command, and on
#     Apple Silicon so does the brew shellenv eval (/opt/homebrew is not on
#     PATH before it); on Linux and WSL2 the apt-get install line, which
#     installs curl, comes before the first curl command;
#   - every NXF_VER= in the page is NEXTFLOW_VERSION from versions.env
#     (a value such as latest fails too).
# Planted copies with one defect each must fail (a check that cannot fail is
# not a check).
#
# Run: bash tests/test_install_docs_order.sh
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
README="${REPO}/README.md"
GUIDE="${REPO}/docs/getting-started.md"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

NXF_PIN=$(sed -n 's/^NEXTFLOW_VERSION="\([0-9.]*\)".*/\1/p' "${REPO}/versions.env")
[ -n "$NXF_PIN" ] || { echo "FAIL: no NEXTFLOW_VERSION in versions.env"; exit 1; }

SECTIONS=("Linux" "Intel Mac" "Apple Silicon Mac" "Windows with WSL2")
VERIFY_CMDS=("docker info" "bash --version" "java -version" "nextflow -version")
# The newest release tag: only v<digit> tags, so fixture-* and e2e-data-* lose.
TAG_PICK="git tag --list 'v[0-9]*' --sort=-v:refname | head -n 1"
# Lines inside a section, as extended regexes.
JAVA_RE='^[[:space:]]*(sudo apt-get install .*openjdk-17|brew install .*temurin@17)'
NXF_RE='NXF_VER='
BREW_RE='^[[:space:]]*brew '
SHELLENV_RE='^[[:space:]]*eval "\$\(/opt/homebrew/bin/brew shellenv\)"'
APT_RE='^[[:space:]]*sudo apt-get install '
CURL_RE='^[[:space:]]*curl '

# section FILE TITLE: the lines under "### TITLE" up to the next heading of
# level 1 to 3 outside a code block.
section() {
  awk -v h="### $2" '
    /^[ \t]*```/ { fence = !fence }
    !fence && $0 == h { on = 1; next }
    on && !fence && /^#{1,3} / { exit }
    on { print }
  ' "$1"
}

# last_block: the contents of the last fenced code block on stdin.
last_block() {
  awk '
    /^[ \t]*```/ { if (fence) { fence = 0 } else { fence = 1; n++; delete b; k = 0 }; next }
    fence { b[++k] = $0 }
    END { for (i = 1; i <= k; i++) print b[i] }
  '
}

# first_block: the contents of the first fenced code block on stdin.
first_block() {
  awk '/^[ \t]*```/ { if (fence) exit; fence = 1; next } fence { print }'
}

# first_line REGEX: the number of the first line on stdin matching the
# extended REGEX, or nothing.
first_line() {
  grep -n -E -- "$1" | head -n 1 | cut -d: -f1 || true
}

# check_checkout WHERE BLOCK: BLOCK clones the repository and checks out the
# newest release tag with the release filter. Prints each problem.
check_checkout() {
  local where=$1 block=$2 rc=0
  if ! grep -q 'git clone https://github.com/GeiserX/OwnGenome' <<< "$block"; then
    echo "${where} does not clone the repository"; rc=1
  fi
  # shellcheck disable=SC2016 # the literal command the docs show
  if ! grep -q 'git checkout "$(git tag' <<< "$block"; then
    echo "${where} does not check out the newest release tag"; rc=1
  elif ! grep -qF "git checkout \"\$(${TAG_PICK})\"" <<< "$block"; then
    echo "${where} picks the tag without the release filter (want: ${TAG_PICK})"; rc=1
  fi
  return "$rc"
}

# check_readme FILE: prints each problem; returns 1 when there is one.
check_readme() {
  local f=$1 rc=0 first_fence begin before block word
  first_fence=$(grep -n '^```' "$f" | head -n 1 | cut -d: -f1)
  begin=$(grep -n -i 'before you begin' "$f" | head -n 1 | cut -d: -f1)
  if [ -z "$first_fence" ]; then
    echo "no fenced code block"; return 1
  fi
  if [ -z "$begin" ] || [ "$begin" -ge "$first_fence" ]; then
    echo "no \"Before you begin\" list before the first code block (line ${first_fence})"
    rc=1
  else
    before=$(sed -n "${begin},${first_fence}p" "$f")
    for word in Docker bash Java Nextflow RAM; do
      grep -q "$word" <<< "$before" || { echo "\"Before you begin\" does not name ${word}"; rc=1; }
    done
    if ! grep 'RAM' <<< "$before" | grep -q 'swap'; then
      echo "the RAM line does not say a 16 GB machine needs swap"; rc=1
    fi
    if ! grep 'RAM' <<< "$before" | grep -q '25-prs/#memory'; then
      echo "the RAM line does not link step 25's memory section"; rc=1
    fi
  fi
  block=$(awk -v s="$first_fence" 'NR > s && /^[ \t]*```/ { exit } NR > s { print }' "$f")
  check_checkout "the first code block" "$block" || rc=1
  return "$rc"
}

# check_guide FILE: prints each problem; returns 1 when there is one.
check_guide() {
  local f=$1 rc=0 t body block cmd brew_line hb_line env_line java_line nxf_line apt_line curl_line v
  if ! grep -E '^\| \*\*RAM\*\*' "$f" | grep -q 'swap'; then
    echo "the RAM row does not say a 16 GB machine needs swap"; rc=1
  fi
  if ! grep -E '^\| \*\*RAM\*\*' "$f" | grep -q '25-prs.md#memory'; then
    echo "the RAM row does not link 25-prs.md#memory"; rc=1
  fi
  body=$(section "$f" "Get the code")
  if [ -z "$body" ]; then
    echo "no \"### Get the code\" section"; rc=1
  else
    check_checkout "Get the code: its first code block" "$(first_block <<< "$body")" || rc=1
  fi
  for t in "${SECTIONS[@]}"; do
    body=$(section "$f" "$t")
    if [ -z "$body" ]; then
      echo "no \"### ${t}\" section"; rc=1; continue
    fi
    grep -q 'Docker' <<< "$body" || { echo "${t}: does not install Docker"; rc=1; }
    grep -qE 'openjdk-17|temurin@17' <<< "$body" || { echo "${t}: does not install Java"; rc=1; }
    grep -q "NXF_VER=" <<< "$body" || { echo "${t}: does not install the pinned Nextflow"; rc=1; }
    block=$(last_block <<< "$body")
    for cmd in "${VERIFY_CMDS[@]}"; do
      grep -qF "$cmd" <<< "$block" || { echo "${t}: the last code block does not run '${cmd}'"; rc=1; }
    done
    if ! awk '/^[ \t]*```/ { fence = !fence; seen = 1; after = 0; next } !fence && seen && /Next/ { after = 1 } END { exit !after }' <<< "$body"; then
      echo "${t}: no \"Next\" line after the verification block"; rc=1
    fi
    # Line numbers below count from the top of the section.
    java_line=$(first_line "$JAVA_RE" <<< "$body")
    nxf_line=$(first_line "$NXF_RE" <<< "$body")
    if [ -z "$java_line" ]; then
      echo "${t}: no apt-get openjdk-17 or brew temurin@17 install line"; rc=1
    elif [ -n "$nxf_line" ] && [ "$nxf_line" -lt "$java_line" ]; then
      echo "${t}: installs Nextflow (line ${nxf_line} of the section) before Java (line ${java_line})"; rc=1
    fi
    case "$t" in
      *Mac)
        brew_line=$(first_line "$BREW_RE" <<< "$body")
        hb_line=$(first_line 'raw\.githubusercontent\.com/Homebrew/install' <<< "$body")
        if [ -z "$hb_line" ]; then
          echo "${t}: no Homebrew install line"; rc=1
        elif [ -n "$brew_line" ] && [ "$brew_line" -lt "$hb_line" ]; then
          echo "${t}: runs brew (line ${brew_line} of the section) before installing Homebrew (line ${hb_line})"; rc=1
        fi
        ;;
      Linux|*WSL2)
        apt_line=$(first_line "$APT_RE" <<< "$body")
        curl_line=$(first_line "$CURL_RE" <<< "$body")
        if [ -z "$apt_line" ]; then
          echo "${t}: no apt-get install line"; rc=1
        elif [ -n "$curl_line" ] && [ "$curl_line" -lt "$apt_line" ]; then
          echo "${t}: runs curl (line ${curl_line} of the section) before apt-get installs it (line ${apt_line})"; rc=1
        fi
        ;;
    esac
    if [ "$t" = "Apple Silicon Mac" ]; then
      env_line=$(first_line "$SHELLENV_RE" <<< "$body")
      brew_line=$(first_line "$BREW_RE" <<< "$body")
      if [ -z "$env_line" ]; then
        echo "${t}: never puts /opt/homebrew on the PATH (no eval of brew shellenv)"; rc=1
      elif [ -n "$brew_line" ] && [ "$brew_line" -lt "$env_line" ]; then
        echo "${t}: runs brew (line ${brew_line} of the section) before brew shellenv puts it on the PATH (line ${env_line})"; rc=1
      fi
    fi
  done
  while read -r v; do
    [ "$v" = "$NXF_PIN" ] || { echo "NXF_VER=${v} is not NEXTFLOW_VERSION ${NXF_PIN}"; rc=1; }
  done < <(grep -oE 'NXF_VER=[^[:space:]]*' "$f" | cut -d= -f2)
  return "$rc"
}

# expect_fail NAME CHECK FILE WANT: a planted copy must fail CHECK, and the
# problems it prints must include WANT (so it fails for the planted defect).
expect_fail() {
  local name=$1 check=$2 file=$3 want=$4 out
  if out=$("$check" "$file"); then
    fail "control: ${name} passed"
  elif ! grep -qF -- "$want" <<< "$out"; then
    fail "control: ${name} failed, but not with \"${want}\": ${out}"
  else
    pass "control: ${name} fails: $(grep -F -- "$want" <<< "$out" | head -n 1)"
  fi
}

# move_after FILE TITLE LINE_RE ANCHOR_RE: FILE with, inside "### TITLE", the
# lines matching LINE_RE taken out and put back after the first later line
# matching ANCHOR_RE (a planted step-order defect).
move_after() {
  awk -v h="### $2" -v pat="$3" -v anc="$4" '
    $0 == h { on = 1; print; next }
    /^[ \t]*```/ { fence = !fence }
    on && !fence && /^#{1,3} / { on = 0 }
    on && !done && $0 ~ pat { held = held $0 "\n"; next }
    on && !done && held != "" && $0 ~ anc { print; printf "%s", held; done = 1; next }
    { print }
  ' "$1"
}

# The real files.
if out=$(check_readme "$README"); then
  pass "README.md lists the requirements before the first command and checks out a release"
else
  fail "README.md:"$'\n'"${out}"
fi
if out=$(check_guide "$GUIDE"); then
  pass "docs/getting-started.md has the four ordered install sections"
else
  fail "docs/getting-started.md:"$'\n'"${out}"
fi

# Controls on README.md.
# shellcheck disable=SC2016 # backticks of a Markdown fence, not a command
{ printf '```bash\necho first\n```\n\n'; cat "$README"; } > "${WORK}/readme-early-block.md"
expect_fail "a code block above \"Before you begin\"" check_readme "${WORK}/readme-early-block.md" 'before the first code block'
# shellcheck disable=SC2016 # sed matches the literal $( of the README line
sed '/git checkout "\$(git tag/d' "$README" > "${WORK}/readme-main.md"
expect_fail "a clone that stays on main" check_readme "${WORK}/readme-main.md" 'newest release tag'
sed 's/--list .v\[0-9\]\*. //' "$README" > "${WORK}/readme-nofilter.md"
expect_fail "a release checkout without the v[0-9]* filter" check_readme "${WORK}/readme-nofilter.md" 'the first code block picks the tag without the release filter'
sed '/RAM/s/swap/memory/g' "$README" > "${WORK}/readme-noswap.md"
expect_fail "a RAM line with no swap" check_readme "${WORK}/readme-noswap.md" 'needs swap'

# Controls on docs/getting-started.md.
sed '/^| \*\*RAM\*\*/s/swap/memory/g' "$GUIDE" > "${WORK}/guide-noswap.md"
expect_fail "a RAM row with no swap" check_guide "${WORK}/guide-noswap.md" 'RAM row does not say'
sed 's/^### Intel Mac$/### Intel Macs/' "$GUIDE" > "${WORK}/guide-nosection.md"
expect_fail "a missing Intel Mac section" check_guide "${WORK}/guide-nosection.md" 'no "### Intel Mac" section'
awk '$0 == "### Apple Silicon Mac" { on = 1 } on && /^nextflow -version$/ { on = 0; next } { print }' "$GUIDE" > "${WORK}/guide-noverify.md"
expect_fail "a verification block without nextflow -version" check_guide "${WORK}/guide-noverify.md" "Apple Silicon Mac: the last code block does not run 'nextflow -version'"
awk '$0 == "### Intel Mac" { on = 1 } on && /^```bash$/ && !done { print; print "brew install wget"; done = 1; next } { print }' "$GUIDE" > "${WORK}/guide-brew-first.md"
expect_fail "a Mac section that runs brew before installing Homebrew" check_guide "${WORK}/guide-brew-first.md" 'Intel Mac: runs brew'
sed 's/--list .v\[0-9\]\*. //' "$GUIDE" > "${WORK}/guide-nofilter.md"
expect_fail "a Get the code checkout without the v[0-9]* filter" check_guide "${WORK}/guide-nofilter.md" 'Get the code: its first code block picks the tag without the release filter'
move_after "$GUIDE" "Intel Mac" 'temurin@17' 'NXF_VER=' > "${WORK}/guide-java-late.md"
expect_fail "a section that installs Nextflow before Java" check_guide "${WORK}/guide-java-late.md" 'Intel Mac: installs Nextflow'
move_after "$GUIDE" "Apple Silicon Mac" 'brew shellenv' 'brew install --cask' > "${WORK}/guide-shellenv-late.md"
expect_fail "an Apple Silicon section that runs brew before brew shellenv" check_guide "${WORK}/guide-shellenv-late.md" 'before brew shellenv puts it on the PATH'
move_after "$GUIDE" "Linux" 'apt-get install' 'get.docker.com' > "${WORK}/guide-curl-first.md"
expect_fail "a Linux section that runs curl before apt-get installs it" check_guide "${WORK}/guide-curl-first.md" 'Linux: runs curl'
awk -v pin="NXF_VER=${NXF_PIN}" '!done && index($0, pin) { sub(pin, "NXF_VER=1.0.0"); done = 1 } { print }' "$GUIDE" > "${WORK}/guide-oldnxf.md"
expect_fail "an NXF_VER that is not the pin" check_guide "${WORK}/guide-oldnxf.md" 'NXF_VER=1.0.0 is not'
awk -v pin="NXF_VER=${NXF_PIN}" '!done && index($0, pin) { sub(pin, "NXF_VER=latest"); done = 1 } { print }' "$GUIDE" > "${WORK}/guide-latestnxf.md"
expect_fail "an unpinned NXF_VER=latest" check_guide "${WORK}/guide-latestnxf.md" 'NXF_VER=latest is not'

if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All install-docs checks passed."
