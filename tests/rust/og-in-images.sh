#!/usr/bin/env bash
# og-in-images.sh: the og binaries in bin/.og run by bare name inside the
# images the pipeline's tasks use, as a Nextflow task calls them: bin/ mounted
# at its own path and first on PATH, the command run by sh. Each image runs on
# every platform it publishes (linux/amd64, linux/arm64; arm64 under QEMU on
# an amd64 host), and og must report the build for that platform's uname -m.
# Then the shim with no binary beside it must stop with exit 127 and the
# setup command.
#
# Usage: tests/rust/og-in-images.sh [IMAGE_VAR ...]
#   default: PYTHON_IMAGE BCFTOOLS_IMAGE MANTA_IMAGE PARASCOPY_IMAGE
# Needs bin/.og built (scripts/setup.sh --og). Exit 1 when any run fails.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
# shellcheck source=../../versions.env
. "${ROOT}/versions.env"
[ $# -gt 0 ] || set -- PYTHON_IMAGE BCFTOOLS_IMAGE MANTA_IMAGE PARASCOPY_IMAGE
BIN="${ROOT}/bin"
for a in x86_64 aarch64; do
  [ -x "${BIN}/.og/og-linux-${a}" ] || { echo "ERROR: ${BIN}/.og/og-linux-${a} is missing: run scripts/setup.sh --og" >&2; exit 1; }
done
fail=0 runs=0
run_og() {  # run_og IMAGE PLATFORM BINDIR: og version, as a task would call it
  # -q: a first pull prints its progress on the output compared below.
  docker run -q --rm --pull missing --platform "$2" -v "${3}:${3}:ro" -e "PATH=${3}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    --entrypoint sh "$1" -c 'echo "$(uname -m) $(og version)"' 2>&1
}
printf '%-18s %-12s %-8s %s\n' image platform result "uname -m and og version"
for var in "$@"; do
  img=${!var:?$var is not set in versions.env}
  for plat in linux/amd64 linux/arm64; do
    out=$(run_og "$img" "$plat" "$BIN") || true
    if [[ "$out" == *"does not provide the specified platform"* || "$out" == *"no matching manifest"* ]]; then
      printf '%-18s %-12s %-8s %s\n' "$var" "$plat" skipped "the image has no ${plat} build"
      continue
    fi
    runs=$((runs + 1))
    want=$([ "$plat" = linux/amd64 ] && echo x86_64 || echo aarch64)
    if [[ "$out" == "${want} og rust tree "*" (${want}-unknown-linux-musl)" ]]; then r=ok; else r=FAIL; fail=1; fi
    printf '%-18s %-12s %-8s %s\n' "$var" "$plat" "$r" "$out"
  done
done
[ "$runs" -gt 0 ] || { echo "ERROR: no image ran"; exit 1; }
# The shim alone: no bin/.og beside it.
# Under rust/target (ignored by git), not TMPDIR: Docker on macOS shares only some folders.
mkdir -p "${ROOT}/rust/target"
empty=$(mktemp -d "${ROOT}/rust/target/og-shim.XXXXXX")
trap 'rm -rf "$empty"' EXIT
cp "${BIN}/og" "${empty}/og"
rc=0
out=$(docker run --rm -v "${empty}:${empty}:ro" -e "PATH=${empty}:/usr/local/bin:/usr/bin:/bin" --entrypoint sh "${!1}" -c 'og version' 2>&1) || rc=$?
if [ "$rc" -eq 127 ] && [[ "$out" == *"og is not built for "*"run scripts/setup.sh --og"* ]]; then
  echo "ok: the shim without a binary exits 127: ${out}"
else
  echo "FAIL: the shim without a binary exited ${rc}: ${out}"; fail=1
fi
[ "$fail" -eq 0 ] && echo "PASS: og runs by bare name in every image and platform above" || echo "FAIL: see above"
exit "$fail"
