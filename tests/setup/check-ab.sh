#!/usr/bin/env bash
# tests/setup/check-ab.sh: run validate-setup.sh's checks twice for every
# call a test makes, with bash (side A, OG_SETUP=bash) and with og (side B,
# og setup check, as validate-setup.sh runs it by default), and compare what
# each side printed on stdout and stderr, its exit code, and the calls it
# made to the fake engine (FAKE_DOCKER_LOG).
#
# Usage: tests/setup/check-ab.sh [--control] -- COMMAND [ARGS...]
#
# COMMAND runs from a copy of this checkout (the tracked and untracked files,
# and bin/.og) whose scripts/validate-setup.sh is a shim: the shim runs the
# real script once per side and then gives the caller side B's output and
# exit code, so COMMAND's own checks read og's output, and keeps side A's
# engine calls in FAKE_DOCKER_LOG (equal to side B's when the two agree).
# For example:
#   tests/setup/check-ab.sh -- scripts/ci/run-fake-docker-suite.sh base-validate
#   tests/setup/check-ab.sh -- bash tests/test_og_check.sh
#
# Side A is validate-setup.sh with one change made in the copy only: it asks
# the engine's memory of CONTAINER_ENGINE, with podman's {{.Host.MemTotal}}
# when {{.MemTotal}} gives no number, as og does (validate-setup.sh asked a
# literal docker). Side B runs with OG_SETUP unset, the default: og when it
# can run og setup check, else validate-setup.sh's own checks; the report
# counts the calls og answered and those it left to bash. A call made with
# OG_SETUP=bash (a test of the bash checks themselves) runs side A only and
# is not counted.
#
# --control adds one byte to side B's stdout in every call before the
# compare: every call must then be reported as differing, or the compare
# could not fail.
#
# Exit 0 when COMMAND exits 0, made at least one call, and every call
# agreed (with --control: every call differed); 1 otherwise; 2 for a wrong
# argument. The copy is kept in CHECK_AB_DIR (default a new temporary
# folder), with each call's outputs in .check-ab/run.N.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
control=false
while [ $# -gt 0 ]; do
  case "$1" in
    --control) control=true; shift ;;
    --) shift; break ;;
    *) echo "usage: tests/setup/check-ab.sh [--control] -- COMMAND [ARGS...]" >&2; exit 2 ;;
  esac
done
[ $# -gt 0 ] || { echo "usage: tests/setup/check-ab.sh [--control] -- COMMAND [ARGS...]" >&2; exit 2; }

dir=${CHECK_AB_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/check-ab.XXXXXX")}
mkdir -p "$dir"
dir=$(cd "$dir" && pwd)
[ -z "$(ls -A "$dir")" ] || { echo "ERROR: ${dir} is not empty" >&2; exit 2; }
copy="${dir}/repo"
ab="${dir}/.check-ab"
mkdir -p "$copy" "$ab"
( cd "$ROOT" && git ls-files -co --exclude-standard -z | while IFS= read -r -d '' f; do [ -e "$f" ] && printf '%s\0' "$f"; done \
  | tar --null -T - -cf - ) | tar -C "$copy" -xf -
if [ -d "${ROOT}/bin/.og" ]; then mkdir -p "${copy}/bin/.og" && cp -p "${ROOT}"/bin/.og/og-* "${copy}/bin/.og/" 2>/dev/null || true; fi
$control && : > "${ab}/control"
# A checkout of its own, so og.sh can stamp and check rust/'s tree.
( cd "$copy" && git init -q && git add -A -f \
  && git -c user.name=check-ab -c user.email=check-ab@example.invalid commit -qm copy )

real="${copy}/scripts/validate-setup.real.sh"
mv "${copy}/scripts/validate-setup.sh" "$real"
# Side A asks CONTAINER_ENGINE for its memory, as og does.
old="DOCKER_MEM_BYTES=\$(docker info --format '{{.MemTotal}}' 2>/dev/null || echo 0)"
grep -qF "$old" "$real" || { echo "ERROR: validate-setup.sh no longer reads the engine memory as check-ab.sh expects" >&2; exit 2; }
new="DOCKER_MEM_BYTES=\$(\"\${CONTAINER_ENGINE:-docker}\" info --format '{{.MemTotal}}' 2>/dev/null || echo 0); if ! [ \"\$DOCKER_MEM_BYTES\" -gt 0 ] 2>/dev/null; then ab_mem=\$(\"\${CONTAINER_ENGINE:-docker}\" info --format '{{.Host.MemTotal}}' 2>/dev/null || echo 0); if [ \"\$ab_mem\" -gt 0 ] 2>/dev/null; then DOCKER_MEM_BYTES=\$ab_mem; fi; fi"
OLD=$old NEW=$new awk 'index($0, ENVIRON["OLD"]) { i = index($0, ENVIRON["OLD"]); $0 = substr($0, 1, i - 1) ENVIRON["NEW"] substr($0, i + length(ENVIRON["OLD"])) } 1' \
  "$real" > "${real}.tmp" && cat "${real}.tmp" > "$real" && rm "${real}.tmp"

cat > "${copy}/scripts/validate-setup.sh" <<SHIM
#!/usr/bin/env bash
# check-ab.sh's shim: side A, side B, compare, then side B's output.
set -uo pipefail
here=\$(cd "\$(dirname "\$0")" && pwd)
# A caller that asks for the bash checks gets them, with no compare.
[ "\${OG_SETUP:-}" = bash ] && exec "\${here}/validate-setup.real.sh" "\$@"
ab='${ab}'
n=1; while ! mkdir "\${ab}/run.\${n}" 2>/dev/null; do n=\$((n + 1)); done
d="\${ab}/run.\${n}"
printf '%q ' "\$@" > "\${d}/args"; pwd > "\${d}/cwd"
log=\${FAKE_DOCKER_LOG:-}
size() { if [ -n "\$log" ] && [ -f "\$log" ]; then wc -c < "\$log" | tr -d ' '; else echo 0; fi; }
pre=\$(size)
OG_SETUP=bash "\${here}/validate-setup.real.sh" "\$@" > "\${d}/a.out" 2> "\${d}/a.err" < /dev/null; echo \$? > "\${d}/a.rc"
mid=\$(size)
og=\${OG:-\${here}/../bin/og}
if "\$og" setup check --help > /dev/null 2>&1; then echo og > "\${d}/b.by"; else echo bash > "\${d}/b.by"; fi
(unset OG_SETUP; "\${here}/validate-setup.real.sh" "\$@") > "\${d}/b.out" 2> "\${d}/b.err" < /dev/null; echo \$? > "\${d}/b.rc"
post=\$(size)
if [ -n "\$log" ] && [ -f "\$log" ]; then
  tail -c +\$((pre + 1)) "\$log" | head -c \$((mid - pre)) > "\${d}/a.log"
  tail -c +\$((mid + 1)) "\$log" > "\${d}/b.log"
  head -c "\$mid" "\$log" > "\${d}/log.keep" && cat "\${d}/log.keep" > "\$log" && rm "\${d}/log.keep"
else
  : > "\${d}/a.log"; : > "\${d}/b.log"
fi
[ -f "\${ab}/control" ] && printf 'X' >> "\${d}/b.out"
same=true
for k in out err rc log; do cmp -s "\${d}/a.\${k}" "\${d}/b.\${k}" || { same=false; echo "\$k" >> "\${d}/differs"; }; done
\$same && : > "\${d}/same"
cat "\${d}/b.out"; cat "\${d}/b.err" >&2
exit "\$(cat "\${d}/b.rc")"
SHIM
chmod +x "${copy}/scripts/validate-setup.sh"

rc=0
(cd "$copy" && "$@") || rc=$?

runs=0 same=0 by_og=0
for d in "${ab}"/run.*; do
  [ -d "$d" ] || continue
  runs=$((runs + 1))
  [ "$(cat "${d}/b.by" 2>/dev/null)" = og ] && by_og=$((by_og + 1))
  if [ -f "${d}/same" ]; then
    same=$((same + 1))
  elif ! $control; then
    echo "=== check-ab $(basename "$d"): validate-setup.sh $(cat "${d}/args")(from $(cat "${d}/cwd")) differs in: $(tr '\n' ' ' < "${d}/differs")"
    for k in $(cat "${d}/differs"); do diff "${d}/a.${k}" "${d}/b.${k}" | head -40 || true; done
  fi
done
echo "check-ab: ${runs} call(s), ${same} the same on both sides, $((runs - same)) different; og answered ${by_og}, bash $((runs - by_og)) (og could not run)."
echo "check-ab: COMMAND exited ${rc}; outputs kept in ${ab}"
if $control; then
  [ "$runs" -gt 0 ] && [ "$same" -eq 0 ] && { echo "check-ab: control: every call was reported as differing: PASS"; exit 0; }
  echo "check-ab: control: a call with a planted difference was reported the same: FAIL"; exit 1
fi
[ "$rc" -eq 0 ] && [ "$runs" -gt 0 ] && [ "$same" -eq "$runs" ] && exit 0
exit 1
