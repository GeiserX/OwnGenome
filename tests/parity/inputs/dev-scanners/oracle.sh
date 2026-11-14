#!/usr/bin/env bash
# Adapt the deleted shell heredocs to the og argv. Run only by parity.
set -euo pipefail
verb=$1 source=$2
shift 2
heredoc() {
  awk -v fn="$1() {" '
    $0 == fn { scope=1 }
    scope && /<<\047(PY|EOF)\047/ { body=1; next }
    body && /^(PY|EOF)$/ { exit }
    body { print }' "$source"
}
inline() {
  awk -v fn="$1() {" '
    $0 == fn { scope=1 }
    scope && /python3 -c \047$/ { body=1; next }
    body && /^\047/ { exit }
    body { print }' "$source"
}
case "$verb" in
  personal)
    root= allow= skip= patterns= cached=
    while [ "$#" -gt 0 ]; do
      case "$1" in
        --root) root=$(cd "$2" && pwd); shift 2 ;;
        --allow) allow=$2; shift 2 ;;
        --skip) skip=$2; shift 2 ;;
        --patterns) patterns=$2; shift 2 ;;
        --cached) cached=1; shift ;;
        --) shift; break ;;
        *) break ;;
      esac
    done
    exec python3 -c "$(heredoc scan)" "$root" "$allow" "$skip" "$patterns" "$cached" "$@" ;;
  images)
    root= exempt= pairs=
    while [ "$#" -gt 0 ]; do
      case "$1" in
        --root) root=$(cd "$2" && pwd); shift 2 ;;
        --exempt-files) exempt=$2; shift 2 ;;
        --exempt-pairs) pairs=$2; shift 2 ;;
        *) exit 2 ;;
      esac
    done
    exec python3 -c "$(heredoc scan)" "$root" "$exempt" "$pairs" "$(cat)" ;;
  coupled) shift; exec python3 -c "$(heredoc coupled)" "$1" ;;
  check-run) shift; exec python3 -c "$(heredoc check_versions)" "$@" ;;
  expand) shift; exec python3 -c "$(inline expand_vars)" "$1" ;;
  trait) shift; export ID=$1; exec python3 -c "$(inline trait_of)" ;;
  label) exec python3 -c "$(heredoc label_matches)" "$@" ;;
  metric)
    while [ "$#" -gt 0 ]; do
      case "$1" in
        --id) export ID=$2 ;;
        --ppm) export PPM=$2 ;;
        --ancestry) export ANC=$2 ;;
        --effect) export EFF=$2 ;;
        *) exit 2 ;;
      esac
      shift 2
    done
    exec python3 -c "$(inline metric_matches)" ;;
  *) exit 2 ;;
esac
