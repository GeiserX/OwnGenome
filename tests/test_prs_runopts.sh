#!/usr/bin/env bash
# The outer run's docker flags reach pgsc_calc's own containers.
#
# PRS (modules/local/prs/main.nf) runs pgsc_calc, a Nextflow pipeline of its
# own, with a config it writes, images.config. Its docker.runOptions line
# must start with -u UID:GID and --network none (and the sex folder's mount),
# then carry the outer run's docker.runOptions and the PRS process's
# containerOptions, less their user and network flags.
#
#   bash tests/test_prs_runopts.sh [MODULE]
#       No Docker, no Nextflow (guard.yml runs it on every pull request).
#       Runs the task script's "outer docker flags" block, with Groovy's
#       escapes undone, on several outer values and reads the line back from
#       the images.config it wrote. MODULE defaults to the PRS module.
#       Without the block (the module before it existed) the test fails.
#
# The Groovy side drops each user flag with its whole value first, by a
# pattern between markers in the module; the test applies that pattern
# before the block. What the Groovy side hands the block (the run's
# docker.runOptions and task.containerOptions) is checked on real runs by the
# e2e case tests/e2e/prs-runopts-and-trace-1.sh.
#
# It also checks the docker profile itself (nextflow.config), whose
# docker.runOptions every task gets: expanded by bash, as a task's
# .command.run expands it, its user is the caller's $(id -u):$(id -g), as
# with the bash steps' run_in; and the block hands it on with no word of
# that user left. The old --user root:root must fail the first check.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
MODULE=${1:-${REPO}/modules/local/prs/main.nf}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0

# The block between the markers, as bash sees it: Groovy turns \\ into \ and
# \$ into $; any other escape in the block is refused, so it cannot be misread.
# shellcheck disable=SC2016  # the dollar signs are awk's and the block's
if ! awk '
  { text = text $0 "\n" }
  END {
    begin = "# >>> outer docker flags"; end = "# <<< outer docker flags"
    b = index(text, begin)
    if (!b) { print "no '\''" begin "'\'' ... '\''" end "'\'' block in " FILENAME > "/dev/stderr"; exit 1 }
    rest = substr(text, b + length(begin))
    e = index(rest, end)
    if (!e) { print "no '\''" begin "'\'' ... '\''" end "'\'' block in " FILENAME > "/dev/stderr"; exit 1 }
    block = substr(rest, 1, e - 1)
    nl = index(block, "\n")
    block = nl ? substr(block, nl + 1) : ""
    out = ""
    for (i = 1; i <= length(block); i++) {
      c = substr(block, i, 1)
      if (c == "\\") {
        nxt = substr(block, i + 1, 1)
        if (nxt != "\\" && nxt != "$") { print "unexpected Groovy escape \\" nxt " in the block" > "/dev/stderr"; exit 1 }
        out = out nxt; i++
      } else out = out c
    }
    print out
  }' "$MODULE" > "${WORK}/block.sh"
then
  echo "FAIL: the outer docker flags block of ${MODULE} could not be read"
  exit 1
fi

# The Groovy side first drops each user flag with its whole value: a value
# can hold spaces ($(id -u):$(id -g)), and the block splits the options into
# words, so the block alone drops only the value's first word. The pattern is
# the single-quoted Groovy string between the module's "user flag" markers, a
# Java regex. perl applies it here as the module's replaceAll(user_flag, ' ')
# does; the pattern uses only syntax the two read the same way. A module
# without the markers hands the options on as they are, as before it dropped
# them, so the cases below fail there.
# shellcheck disable=SC2016  # the dollar signs are awk's
if ! USER_FLAG=$(awk '
  /\/\/ >>> user flag/ { on = 1; next }
  /\/\/ <<< user flag/ { on = 0 }
  on && /^ *def user_flag = '\''.*'\'' *$/ {
    s = $0; sub(/^ *def user_flag = '\''/, "", s); sub(/'\'' *$/, "", s)
    out = ""
    for (i = 1; i <= length(s); i++) {
      c = substr(s, i, 1)
      if (c == "\\") {
        nxt = substr(s, i + 1, 1)
        if (nxt != "\\") { print "unexpected Groovy escape \\" nxt " in user_flag" > "/dev/stderr"; exit 1 }
        out = out nxt; i++
      } else out = out c
    }
    print out; n++
  }
  END { if (n > 1) { print "more than one user_flag in " FILENAME > "/dev/stderr"; exit 1 } }' "$MODULE"); then
  echo "FAIL: the user flag pattern of ${MODULE} could not be read"
  exit 1
fi
echo "user flag pattern: ${USER_FLAG:-none, the options pass as they are}"
# groovy_side OUTER: OUTER after the module's replaceAll(user_flag, ' ').
groovy_side() {
  if [ -z "$USER_FLAG" ]; then printf '%s' "$1"; return; fi
  # shellcheck disable=SC2016  # perl's own variable
  printf '%s' "$1" | USER_FLAG=$USER_FLAG perl -0777 -pe 's/$ENV{USER_FLAG}/ /g'
}

# line OUTER MOUNT: the docker.runOptions line the block writes for those
# values, after the process block the task writes first (a file a glob would
# find). OUTER goes through the Groovy side first, as in the module.
line() {
  local d outer
  d=$(mktemp -d "${WORK}/case.XXXXXX")
  outer=$(groovy_side "$1")
  printf 'process {\n}\n' > "${d}/images.config"
  (cd "$d" && OUTER_DOCKER=$outer MOUNT=$2 bash -euo pipefail "${WORK}/block.sh" > /dev/null) || { echo "the block failed"; return; }
  grep '^docker\.runOptions' "${d}/images.config" || echo "no docker.runOptions line"
}

FIXED="-u $(id -u):$(id -g) --network none"
# expect DESCRIPTION OUTER MOUNT WANT
expect() {
  local got
  got=$(line "$2" "$3")
  if [ "$got" = "docker.runOptions = '$4'" ]; then
    echo "ok:   $1"
  else
    echo "FAIL: $1"
    echo "      outer: $2"
    echo "      got:   ${got}"
    echo "      want:  docker.runOptions = '$4'"
    FAILS=$((FAILS + 1))
  fi
}

expect "no outer flags: the line of before" "" "" "$FIXED"
expect "a config's --user root:root and --network none are dropped" \
  "--user root:root --network none" "" "$FIXED"
expect "a low-priority runOptions reaches pgsc_calc after -u and --network none" \
  "--user root:root --cpu-shares 256" "" "${FIXED} --cpu-shares 256"
expect "every user and network form is dropped, the rest kept in order" \
  "-u 0:0 --cpu-shares 256 --user=1:1 --net host -u0 --network=bridge --memory-swappiness 0 --network none --label a=b" "" \
  "${FIXED} --cpu-shares 256 --memory-swappiness 0 --label a=b"
expect "the sex folder's mount stays before the outer flags" \
  "--cpu-shares 256" ' -v "/w/sex:/w/sex:ro"' "${FIXED}"' -v "/w/sex:/w/sex:ro" --cpu-shares 256'
expect "a glob is passed as it is, not expanded to the files of the task folder" \
  "--cpu-shares 256 *" "" "${FIXED} --cpu-shares 256 *"
# A quote and a backslash, escaped for the single-quoted config string.
expect "a quote and a backslash are escaped for the config" \
  "--label 'a=b c' --env X=a\\b" "" "${FIXED} --label \\'a=b c\\' --env X=a\\\\b"
# A user flag whose value holds spaces is dropped whole, not word by word.
# shellcheck disable=SC2016  # the options' own text, expanded by the task
expect "a config's -u \$(id -u):\$(id -g), the common form, is dropped whole" \
  '-u $(id -u):$(id -g) -e HOME=/tmp' "" "${FIXED} -e HOME=/tmp"
# shellcheck disable=SC2016
expect "--user=\$(...), a quoted --user value and a backquoted -u value are dropped whole" \
  '--user=$(id -u):$(id -g) --cpu-shares 256 --user "$(id -u):$(id -g)" -u `id -u`:`id -g` --label a=b' "" \
  "${FIXED} --cpu-shares 256 --label a=b"
# shellcheck disable=SC2016
expect "the space-free -u \${UID}:\${GROUPS} passes clean" \
  '-u ${UID}:${GROUPS} -e HOME=/tmp' "" "${FIXED} -e HOME=/tmp"
expect "--userns is not a user flag and is kept" \
  "--userns host --cpu-shares 256" "" "${FIXED} --userns host --cpu-shares 256"

# --- The docker profile: every task runs as the caller ---------------------------
# runs_as_caller OPTIONS: OPTIONS, split and expanded by bash as a task's
# .command.run does, give the container the caller's ids. The last -u or
# --user wins, as in `docker run`.
runs_as_caller() {
  local user='' a
  local -a w=()
  eval "w=($1)"
  set -- ${w[@]+"${w[@]}"}
  while [ $# -gt 0 ]; do
    a=$1; shift
    case "$a" in
      -u|--user) user=${1-}; [ $# -eq 0 ] || shift ;;
      -u?*) user=${a#-u} ;;
      --user=*) user=${a#--user=} ;;
    esac
  done
  [ "$user" = "$(id -u):$(id -g)" ]
}
# The docker profile's runOptions, the one docker.runOptions line of nextflow.config.
OPTS=$(sed -n "s/^ *docker\.runOptions *= *'\(.*\)'\$/\1/p" "${REPO}/nextflow.config")
if [ "$(printf '%s\n' "$OPTS" | grep -c .)" -ne 1 ]; then
  echo "FAIL: want one docker.runOptions line in nextflow.config, got: ${OPTS:-none}"
  FAILS=$((FAILS + 1))
elif runs_as_caller "$OPTS"; then
  echo "ok:   the docker profile runs every task as the caller ($(id -u):$(id -g)): ${OPTS}"
else
  echo "FAIL: the docker profile does not run tasks as the caller ($(id -u):$(id -g)): ${OPTS}"
  FAILS=$((FAILS + 1))
fi
# Controls: the profile of before, and the user flag dropped, are refused.
for bad in "--user root:root" "-e HOME=/tmp"; do
  if runs_as_caller "$bad"; then
    echo "FAIL: control: '${bad}' was read as running as the caller"
    FAILS=$((FAILS + 1))
  else
    echo "ok:   control: '${bad}' is refused"
  fi
done
# Through the block, pgsc_calc gets the profile's other flags and none of its user.
expect "the docker profile's own runOptions: its user is dropped, the rest kept" "$OPTS" "" "${FIXED} -e HOME=/tmp"

# Control: a block that drops nothing must fail the second case, so the
# comparison above can fail.
# shellcheck disable=SC2016  # $out and $o belong to the block
sed 's/skip=true ;;/out="$out $o" ;;/' "${WORK}/block.sh" > "${WORK}/broken.sh"
mv "${WORK}/block.sh" "${WORK}/good.sh"
mv "${WORK}/broken.sh" "${WORK}/block.sh"
USER_FLAG=''   # the block alone: the Groovy side would drop the user first
if [ "$(line "--user root:root --network none" "")" = "docker.runOptions = '${FIXED}'" ]; then
  echo "FAIL: control: a block that keeps --user root:root was not caught"
  FAILS=$((FAILS + 1))
else
  echo "ok:   control: a block that keeps --user root:root is caught"
fi

echo "${FAILS} failure(s)"
[ "$FAILS" -eq 0 ]
