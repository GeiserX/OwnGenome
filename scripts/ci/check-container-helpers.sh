#!/usr/bin/env bash
# check-container-helpers.sh: every helper binary (bcftools, bgzip, tabix,
# samtools) that a script or Nextflow module calls must exist in the image the
# call runs in.
#
# The container smoke test (tests/smoke/commands.tsv) runs each image's own
# tool, not the helpers a step calls beside it, so it cannot see a module that
# calls bcftools inside the delly image, or a script that runs bgzip inside
# staphb/bcftools, which ships only bcftools.
#
# How the table is built (no hand-kept list):
#   scripts/*.sh   every `docker run ... IMAGE cmd` command, and every call of
#                  the wrapper from scripts/lib/common.sh,
#                  `run_in [--net] [--root] [--rw DIR] ... IMAGE cmd`; the
#                  in-image commands are `cmd`, or, for `bash -c '...'`, the
#                  first word of every command in the quoted body. A pipe
#                  after the closing quote runs on the host and is not counted.
#                  A tree whose scripts yield no helper call at all fails: the
#                  parser no longer sees how the scripts start containers.
#   modules/local/*/main.nf
#                  every process: its `container` (or its withName selector
#                  in conf/containers.config) and the first word of every
#                  command in the triple-quoted strings of its script:,
#                  shell: and exec: blocks (never stub:). Nextflow runs an
#                  exec: block on the host, so a helper call there has no
#                  image and fails.
# Only calls to the helper binaries above are kept, by bare name: a call by
# path ("$LIBEXEC/bgzip", Manta's own copy) is not a PATH lookup, so
# `command -v` would test another binary. Images are named by their
# versions.env variable when one matches; versions.env is read, never
# executed, and a line that is not NAME="value" stops the check (exit 2).
#
# The table and rules 1 and 3 are og dev helpers
# (rust/crates/og-dev/src/helpers.rs), run through bin/og or the binary OG
# names; this script keeps the docker probe, the orad download and the
# self-test.
#
# Three rules:
#   1. (static) bgzip and tabix are never called in BCFTOOLS_IMAGE or
#      SAMTOOLS_IMAGE: use `bcftools view -Oz -o` and `bcftools index -t`.
#   2. (docker) each image is pulled and `command -v` checks each binary.
#   3. (--orad) orad, Illumina's ORA decompressor, runs on the host, not in an
#      image. The pinned orad release is downloaded, checked against its
#      sha256 and run with --help; every option a script passes to "$ORAD"
#      must be in that list. Step 01 once passed --output-directory, which
#      orad does not have.
#
# Scope: only these four helpers are checked. Probing every in-image command
# (python3, awk, Rscript...) would mean pulling every image the pipeline
# uses, tens of GB, on each pull request.
#
# Usage:
#   scripts/ci/check-container-helpers.sh            both rules (needs docker)
#   scripts/ci/check-container-helpers.sh --static   rule 1 only
#   scripts/ci/check-container-helpers.sh --list     print the derived table
#   scripts/ci/check-container-helpers.sh --orad     rule 3 only (needs curl;
#                                                    --help-file FILE reads
#                                                    the option list from FILE)
#   scripts/ci/check-container-helpers.sh --root DIR check another tree
#   scripts/ci/check-container-helpers.sh --self-test  prove rules 1 and 3
#                                                      can fail, and that a
#                                                      tree whose calls are
#                                                      not recognised fails too
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

if [ "${1:-}" = "--self-test" ]; then
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' EXIT
  for t in bad good; do
    mkdir -p "${tmp}/${t}/scripts" "${tmp}/${t}/modules/local/x"
    printf '%s\n' 'BCFTOOLS_IMAGE="example/bcftools:1.0"' 'SAMTOOLS_IMAGE="example/samtools:1.0"' \
      'TOOL_IMAGE="example/tool:1.0"' > "${tmp}/${t}/versions.env"
  done
  # Planted calls, each behind a different prefix the parser must look past.
  # shellcheck disable=SC2016  # the dollar signs are the test input
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'if ! docker run --rm -v "$G:/genome" "${BCFTOOLS_IMAGE}" bash -c "bcftools view x | bgzip -c > y.gz && tabix y.gz"; then' \
    '  exit 1' \
    'fi' \
    'OUT=$(time docker run --rm "${BCFTOOLS_IMAGE}" tabix z.gz)' \
    'VAR=x docker run --rm "${SAMTOOLS_IMAGE}" bgzip w' \
    'docker run --rm "${TOOL_IMAGE}" bgzip allowed.vcf' \
    'docker run --rm "${BCFTOOLS_IMAGE}" bash -euo pipefail -c "bgzip combined.vcf"' \
    'IMG="${BCFTOOLS_IMAGE}"' \
    'docker run --rm "$IMG" tabix indirect.vcf.gz' \
    'ODD=$(pick_image)' \
    'docker run --rm "$ODD" samtools view x.bam' \
    'run_in --net --root --rw "$G/reference" --cpus 2 "${SAMTOOLS_IMAGE}" tabix wrapped.vcf.gz' \
    'COUNT=$(run_in "${BCFTOOLS_IMAGE}" bash -c "bcftools view x | bgzip -c | wc -l")' \
    > "${tmp}/bad/scripts/a.sh"
  # A withName selector overrides the process container directive.
  mkdir -p "${tmp}/bad/conf" "${tmp}/bad/modules/local/y"
  printf '%s\n' "process { withName: 'Y' { container = 'example/bcftools:1.0' } }" \
    > "${tmp}/bad/conf/containers.config"
  printf '%s\n' 'process Y {' "    container 'example/tool:1.0'" '    script:' '    """' \
    '    tabix y.vcf.gz' '    """' '}' > "${tmp}/bad/modules/local/y/main.nf"
  printf '%s\n' \
    'process X {' \
    "    container 'example/bcftools:1.0'" \
    '    script:' \
    '    """' \
    '    bcftools view in.vcf | bgzip -c > out.vcf.gz' \
    '    """' \
    '    stub:' \
    '    """' \
    '    tabix stub_only.vcf.gz' \
    '    """' \
    '}' > "${tmp}/bad/modules/local/x/main.nf"
  # A shell: block runs its commands in the image too; an exec: block runs
  # on the host, so its helper call has no image and must fail.
  mkdir -p "${tmp}/bad/modules/local/z" "${tmp}/bad/modules/local/w"
  printf '%s\n' 'process Z {' "    container 'example/bcftools:1.0'" '    shell:' "    '''" \
    '    bgzip -c !{vcf} > out.vcf.gz' "    '''" '}' > "${tmp}/bad/modules/local/z/main.nf"
  printf '%s\n' 'process W {' "    container 'example/samtools:1.0'" '    exec:' '    """' \
    '    tabix w.vcf.gz' '    """.execute().waitFor()' '}' > "${tmp}/bad/modules/local/w/main.nf"
  # Nextflow strips a script's common indent before bash runs it, so an
  # indented <<- here-document end closes the body and the call after it runs.
  mkdir -p "${tmp}/bad/modules/local/v"
  printf '%s\n' 'process V {' "    container 'example/bcftools:1.0'" '    script:' '    """' \
    '    cat <<-END_VERSIONS > versions.yml' '    "V": 1' '    END_VERSIONS' \
    '    bgzip -c in.vcf > out.vcf.gz' '    """' '}' > "${tmp}/bad/modules/local/v/main.nf"
  # shellcheck disable=SC2016
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'if ! docker run --rm "${BCFTOOLS_IMAGE}" bash -c "bcftools view x -Oz -o y.gz && bcftools index -t y.gz"; then' \
    '  exit 1' \
    'fi' \
    'docker run --rm "${TOOL_IMAGE}" bgzip allowed.vcf' \
    'docker run --rm "${BCFTOOLS_IMAGE}" bash -euo pipefail -c "bcftools index -t ok.vcf.gz"' \
    'run_in --rw "$G/reference" --cpus 2 "${TOOL_IMAGE}" tabix wrapped.vcf.gz' \
    'run_in "${BCFTOOLS_IMAGE}" bcftools index -t wrapped.vcf.gz' \
    'run_in "${BCFTOOLS_IMAGE}" bash -c "\"\$d/bgzip\" -c by_path.vcf"' \
    > "${tmp}/good/scripts/a.sh"
  # A tree whose scripts start containers in a way the parser does not know:
  # the module alone must not make the check pass.
  mkdir -p "${tmp}/blind/scripts" "${tmp}/blind/modules/local/x"
  cp "${tmp}/good/versions.env" "${tmp}/blind/versions.env"
  # shellcheck disable=SC2016
  printf '%s\n' '#!/usr/bin/env bash' 'start_container "${BCFTOOLS_IMAGE}" bcftools view x' \
    > "${tmp}/blind/scripts/a.sh"
  sed 's/| bgzip -c > out.vcf.gz/-Oz -o out.vcf.gz/' "${tmp}/bad/modules/local/x/main.nf" \
    > "${tmp}/good/modules/local/x/main.nf"

  fail=0
  rc=0; out=$("$0" --static --root "${tmp}/bad" 2>&1) || rc=$?
  [ "$rc" -eq 1 ] || { echo "self-test: planted tree exited ${rc}, expected 1"; fail=1; }
  # shellcheck disable=SC2016  # regexes, not expansions
  for want in 'bgzip +in BCFTOOLS_IMAGE +scripts/a\.sh:2$' 'tabix +in BCFTOOLS_IMAGE +scripts/a\.sh:2$' \
              'tabix +in BCFTOOLS_IMAGE +scripts/a\.sh:5$' 'bgzip +in SAMTOOLS_IMAGE +scripts/a\.sh:6$' \
              'bgzip +in BCFTOOLS_IMAGE +modules/local/x/main\.nf:5 \(X\)$' \
              'bgzip +in BCFTOOLS_IMAGE +scripts/a\.sh:8$' 'tabix +in BCFTOOLS_IMAGE +scripts/a\.sh:10$' \
              'scripts/a\.sh:12 \(\$ODD\)' 'tabix +in BCFTOOLS_IMAGE +modules/local/y/main\.nf:5 \(Y\)$' \
              'tabix +in SAMTOOLS_IMAGE +scripts/a\.sh:13$' 'bgzip +in BCFTOOLS_IMAGE +scripts/a\.sh:14$' \
              'bgzip +in BCFTOOLS_IMAGE +modules/local/z/main\.nf:5 \(Z\)$' \
              'modules/local/w/main\.nf:5 \(W\) exec: runs on the host$' \
              'bgzip +in BCFTOOLS_IMAGE +modules/local/v/main\.nf:8 \(V\)$'; do
    grep -qE -- "$want" <<<"$out" || { echo "self-test: planted call not reported: ${want}"; fail=1; }
  done
  for nope in 'TOOL_IMAGE' 'stub_only' ':9 \(X\)' 'IMAGE +modules/local/w/'; do
    if grep -qE -- "$nope" <<<"$out"; then echo "self-test: reported although allowed: ${nope}"; fail=1; fi
  done
  [ "$fail" -eq 0 ] || printf '%s\n' "$out"
  rc=0; out=$("$0" --static --root "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -ne 0 ]; then
    echo "self-test: clean tree exited ${rc}, expected 0"
    printf '%s\n' "$out"
    fail=1
  fi
  sed 's/| bgzip -c > out.vcf.gz/-Oz -o out.vcf.gz/' "${tmp}/bad/modules/local/x/main.nf" \
    > "${tmp}/blind/modules/local/x/main.nf"
  rc=0; out=$("$0" --static --root "${tmp}/blind" 2>&1) || rc=$?
  if [ "$rc" -ne 2 ] || ! grep -q 'no helper call in scripts/' <<<"$out"; then
    echo "self-test: a tree whose script calls are not recognised exited ${rc}, expected 2 with 'no helper call in scripts/'"
    printf '%s\n' "$out"
    fail=1
  fi
  # Rule 3 against a planted orad option list: an option orad does not list
  # fails, the options it lists pass, and a tree with no orad call or an
  # empty help text fails too.
  printf '%s\n' '-P  path, --path  path   : specify where to write the output file' \
    '-h, --help               : print help and exit' \
    '--ora-reference          : set directory path containing the ora reference file' > "${tmp}/orad-help"
  : > "${tmp}/orad-help-empty"
  # shellcheck disable=SC2016,SC1003  # literal $ and line-continuation backslashes are the test input
  printf '%s\n' '#!/usr/bin/env bash' 'if ! command -v "$ORAD" >/dev/null; then exit 1; fi' \
    '"$ORAD" \' '  --ora-reference "$R" \' '  --output-directory "$O" \' '  "$F"' > "${tmp}/bad/scripts/b.sh"
  # shellcheck disable=SC2016
  printf '%s\n' '#!/usr/bin/env bash' 'time "${ORAD}" --ora-reference="$R" -P "$O" --path "$O" "$F"' \
    > "${tmp}/good/scripts/b.sh"
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help" --root "${tmp}/bad" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q 'scripts/b\.sh:3 passes --output-directory' <<<"$out" \
      || ! grep -q 'OK: --ora-reference' <<<"$out"; then
    echo "self-test: an orad option missing from --help exited ${rc}, expected 1 naming scripts/b.sh:3 --output-directory"
    printf '%s\n' "$out"
    fail=1
  fi
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help" --root "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -ne 0 ] || [ "$(grep -c '^OK: ' <<<"$out")" -ne 3 ]; then
    echo "self-test: orad options that --help lists exited ${rc}, expected 0 with 3 OK lines"
    printf '%s\n' "$out"
    fail=1
  fi
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help-empty" --root "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ]; then
    echo "self-test: an empty orad --help exited ${rc}, expected 1"
    printf '%s\n' "$out"
    fail=1
  fi
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help" --root "${tmp}/blind" 2>&1) || rc=$?
  # shellcheck disable=SC2016  # a literal $ORAD in the expected message
  if [ "$rc" -ne 2 ] || ! grep -q 'no "\$ORAD" call' <<<"$out"; then
    echo "self-test: a tree with no orad call exited ${rc}, expected 2"
    printf '%s\n' "$out"
    fail=1
  fi
  if [ "$fail" -eq 0 ]; then
    echo "self-test: bgzip/tabix behind if !, time, VAR=, \$(...), bash -euo -c, an image alias, the run_in wrapper and in the script:, shell: and exec: blocks of modules, also after an indented here-document end, are reported; a clean tree passes; a tree whose script calls are not recognised fails; an orad option missing from --help, an empty --help and a tree with no orad call fail: PASS"
    exit 0
  fi
  echo "self-test: FAIL"
  exit 1
fi

# Everything but the docker probe and the orad download is og dev helpers
# (rust/crates/og-dev/src/helpers.rs). OG names another og binary.
OG=${OG:-"${ROOT}/bin/og"}

# Rule 3's pin. orad is Illumina's own binary, not an image: step 01 calls it
# on the host. Only the first 9 MB of the archive are read: the binary comes
# before the 785 MB reference.
ORAD_VERSION="2.7.0"
ORAD_URL="https://webdata.illumina.com/downloads/software/dragen-decompression/orad.${ORAD_VERSION}.linux.tar.gz"
ORAD_MEMBER="orad.${ORAD_VERSION}.linux/orad"
ORAD_SHA256="2cc4a6bb9d0721c6555e01bb23f3554753166bfdbdd89ee9f4a8639144a1f066"

mode=all root=$ROOT help_file=""
while [ $# -gt 0 ]; do
  case $1 in
    --static) mode=static ;;
    --list) mode=list ;;
    --orad) mode=orad ;;
    --root) root=$(cd "${2:?--root needs a directory}" && pwd); shift ;;
    --help-file) help_file=${2:?--help-file needs a file}; shift ;;
    *) echo "ERROR: unknown argument $1 (see the usage at the top of $0)" >&2; exit 2 ;;
  esac
  shift
done

case $mode in
  list) exec "$OG" dev helpers --list --root "$root" ;;
  static) exec "$OG" dev helpers --static --root "$root" ;;
  orad)
    # No orad call is exit 2 before anything is downloaded.
    "$OG" dev helpers --orad-calls --root "$root"
    if [ -z "$help_file" ]; then
      tmp=$(mktemp -d)
      trap 'rm -rf "$tmp"' EXIT
      stop=-q   # stop reading once the binary is out
      if tar --version 2>/dev/null | grep -q GNU; then stop=--occurrence=1; fi
      # The exit status is not used: tar stops early, so curl ends on a
      # broken pipe. The checksum below is what proves the right binary
      # arrived. --retry covers a transient network error; curl never
      # retries its (23) write error, so the early stop still ends the fetch.
      echo "Fetching ${ORAD_MEMBER} from ${ORAD_URL} (tar stops after it, so a curl (23) write error below is expected)"
      curl -fsSL --retry 3 "$ORAD_URL" | tar -xzf - -C "$tmp" "$stop" "$ORAD_MEMBER" || true
      binary="${tmp}/${ORAD_MEMBER}"
      if [ ! -f "$binary" ]; then
        echo "FAIL: could not fetch ${ORAD_MEMBER} from ${ORAD_URL}"
        exit 1
      fi
      if command -v sha256sum >/dev/null; then sha=$(sha256sum "$binary"); else sha=$(shasum -a 256 "$binary"); fi
      sha=${sha%% *}
      if [ "$sha" != "$ORAD_SHA256" ]; then
        echo "FAIL: orad ${ORAD_VERSION} has sha256 ${sha}, pinned ${ORAD_SHA256}"
        exit 1
      fi
      chmod 755 "$binary"
      help_file="${tmp}/help"
      "$binary" --help > "$help_file" 2>&1 || true
    fi
    exec "$OG" dev helpers --orad --help-file "$help_file" --orad-version "$ORAD_VERSION" --root "$root"
    ;;
esac

# Rule 1, then rule 2: pull each image, run `command -v` for each helper it
# is asked for, remove it.
rc=0
"$OG" dev helpers --static --root "$root" || rc=$?
[ "$rc" -le 1 ] || exit "$rc"
failed=$rc
table=$("$OG" dev helpers --images --root "$root")
n=0 missing_total=0
# shellcheck disable=SC2016  # the probe runs in the image
probe='for b in "$@"; do command -v "$b" >/dev/null 2>&1 && echo "found $b" || echo "missing $b"; done'
while IFS=$'\t' read -r kind name value bins; do
  [ "$kind" = image ] || continue
  n=$((n + 1))
  echo "::group::${name} (${value}): ${bins}"
  if ! err=$(docker pull -q "$value" 2>&1 >/dev/null); then
    echo "::endgroup::"
    echo "FAIL: cannot pull ${value}: ${err}"
    failed=1
    continue
  fi
  prc=0
  # shellcheck disable=SC2086  # one argument per binary
  res=$(docker run --rm --entrypoint sh "$value" -c "$probe" sh $bins 2>&1) || prc=$?
  [ -z "$res" ] || printf '%s\n' "$res"
  echo "::endgroup::"
  if [ "$prc" -ne 0 ] || ! grep -qE '^(found|missing) ' <<<"$res"; then
    echo "FAIL: could not probe ${value} (exit ${prc})"
    failed=1
    continue
  fi
  for b in $bins; do
    if ! grep -qx "found ${b}" <<<"$res"; then
      missing_total=$((missing_total + 1))
      failed=1
      echo "FAIL: ${b} is not in ${name} (${value}), called at:"
      awk -F'\t' -v n="$name" -v v="$value" -v b="$b" \
        '$1 == "call" && $2 == n && $3 == v && $4 == b { print "    " $5 }' <<<"$table"
    fi
  done
  docker rmi -f "$value" >/dev/null 2>&1 || true
done <<<"$table"

if [ "$missing_total" -eq 0 ] && [ "$failed" -eq 0 ]; then
  echo "OK: every helper binary is present in the image that calls it (${n} images)."
fi
exit "$failed"
