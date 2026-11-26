#!/usr/bin/env bash
# RECORD (this commit only). Writes the CPython answers under
# tests/rust/vectors/ from the oracles this commit still has, the way CI
# ran them: tests/rust/pycompat_oracle.py in PYTHON_IMAGE, and the Rust
# tests that call python3 in a Linux container with the CPython 3.9, 3.12
# and 3.15 that uv 0.12.18 installs (guard.yml), 3.15 first on PATH.
# Needs docker and the network (crates and uv).
#
# Usage: tests/rust/vectors/make.sh WORK   (an empty scratch folder)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
OUT="${ROOT}/tests/rust/vectors"
WORK=$(cd "${1:?Usage: $0 WORK}" && pwd)
# shellcheck source=../../../versions.env
. "${ROOT}/versions.env"
RUST_IMAGE="rust:$(sed -n 's/^channel = "\(.*\)"$/\1/p' "${ROOT}/rust/rust-toolchain.toml")-trixie"
mkdir -p "${WORK}/rec" "${OUT}/pycompat" "${OUT}/scanner"

# 1. pycompat: 1000 cases per operation, a fixed seed, and the tree the
#    tree operations walked, listed by GNU find.
docker run --rm --network none -e PYCOMPAT_TREE=/tree "$PYTHON_IMAGE" sh -c \
  'mkdir /tree && python3 -I -c "$1" 1000 20261126 > /tmp/pairs && find /tree -mindepth 1 -printf "%y\t%P\t%l\t%T@\n" > /tmp/tree && find /tree -maxdepth 0 -printf "%y\t\t\t%T@\n" >> /tmp/tree && cat /tmp/pairs && echo @@TREE@@ && cat /tmp/tree' \
  _ "$(cat "${ROOT}/tests/rust/pycompat_oracle.py")" > "${WORK}/oracle.out"
LC_ALL=C sed '/^@@TREE@@$/,$d' "${WORK}/oracle.out" > "${WORK}/pairs.txt"
LC_ALL=C sed '1,/^@@TREE@@$/d' "${WORK}/oracle.out" > "${WORK}/tree.txt"

# 2. The Rust tests, recording what python3 answered, on this commit's
#    files (git archive) in a docker volume. The pinned oracle commits the
#    tests read come from the public repository.
# CARGO_LOCK_DIR, when set, is a mkdir lock held while the build runs.
LOCKDIR=${CARGO_LOCK_DIR:-}
if [ -n "$LOCKDIR" ]; then
  until mkdir "$LOCKDIR" 2>/dev/null; do sleep 30; done
  printf "%s %s\n" "${LOCK_LABEL:-make.sh}" "$(date -u +%FT%TZ)" > "$LOCKDIR/owner"
  trap 'rm -rf "$LOCKDIR"' EXIT
fi
VOL=og-record-vectors
docker volume create "$VOL" > /dev/null
docker run -i --rm -v "${VOL}:/work" "$RUST_IMAGE" sh -c 'rm -rf /work/rec /work/src && mkdir -p /work/rec && cat > /work/pairs.txt' < "${WORK}/pairs.txt"
git -C "$ROOT" archive HEAD | docker run -i --rm -v "${VOL}:/work" \
  -e CARGO_TARGET_DIR=/work/target -e CARGO_HOME=/work/cargo-home -e CARGO_BUILD_JOBS=4 \
  -e CARGO_PROFILE_DEV_DEBUG=0 -e CARGO_PROFILE_TEST_DEBUG=0 -e CARGO_INCREMENTAL=0 \
  -e OG_RECORD_VECTORS=/work/rec -e PYCOMPAT_RAW=/work/pairs.txt -e PYCOMPAT_CHARS=/work/chars.tsv \
  "$RUST_IMAGE" bash -euo pipefail -c '
    mkdir -p /work/src && tar x -C /work/src && cd /work/src
    apt-get update -qq && apt-get install -y -qq jq clang libclang-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev > /dev/null
    git config --global --add safe.directory "*"
    git init -q && git add -A
    git remote add origin https://github.com/GeiserX/OwnGenome.git
    git fetch -q --depth=1 origin 4104a93f2121b457413336574b8563fe85291fd7 630292bd1aef8cb94bbed604e8961a98facdb47d \
      669a8196cf4c377200f9b4af9c235ef85fb774a1 1ec630690e63796b4410d23ec65768200adb8dc0
    curl -LsSf https://astral.sh/uv/0.12.18/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh > /dev/null
    uv python install --no-bin 3.9 3.12 3.15
    interpreters=""
    for version in 3.9 3.12 3.15; do
      python=$(uv python find --managed-python "$version")
      "$python" -I -c "import sys,unicodedata; print(sys.version.split()[0],unicodedata.unidata_version)"
      interpreters="${interpreters:+$interpreters:}$python"
    done
    export PERSONAL_DATA_PYTHONS=$interpreters PERSONAL_DATA_REQUIRE_PYTHONS=3.9:3.12:3.15
    export PATH="$(dirname "$python"):$PATH"
    python3 --version
    cd rust
    step() { "$@" || echo "RECORD STEP FAILED: $*"; }
    step cargo test --locked -p og-dev --lib -- --include-ignored personal_data::tests:: lookaround::python_tests::
    step cargo test --locked -p og-dev --test ci_readers -- --ignored differential_
    step cargo test --locked -p og-runlog --test parity
    step cargo test --locked -p og-live --test ops -- the_card_and_metrics
    step cargo test --locked -p og --no-default-features --bin og -- --ignored make_chars_table
  '
rm -rf "${WORK}/rec" "${WORK}/chars.tsv"
docker run --rm -v "${VOL}:/work" "$RUST_IMAGE" tar c -C /work rec chars.tsv | tar x -C "$WORK"
[ -z "$LOCKDIR" ] || { rm -rf "$LOCKDIR"; trap - EXIT; }

# 3. Assemble: one row per key, sorted, gzipped without a timestamp.
for f in spans model oracle regex_case unicode; do
  LC_ALL=C sort -u "${WORK}/rec/${f}.jsonl" | gzip -9 -n > "${OUT}/scanner/${f}.jsonl.gz"
done
for f in ci_readers runlog live_ops; do
  LC_ALL=C sort -u "${WORK}/rec/${f}.jsonl" | gzip -9 -n > "${OUT}/${f}.jsonl.gz"
done
LC_ALL=C awk 'NR % 2 == 1 { chars = ($0 ~ /^chars\t/) } { print (NR % 2 == 0 && chars) ? "@CHARS@" : $0 }' \
  "${WORK}/pairs.txt" | gzip -9 -n > "${OUT}/pycompat/cases.tsv.gz"
gzip -9 -n < "${WORK}/chars.tsv" > "${OUT}/pycompat/chars.tsv.gz"
gzip -9 -n < "${WORK}/tree.txt" > "${OUT}/pycompat/tree.tsv.gz"
ls -l "${OUT}" "${OUT}/pycompat" "${OUT}/scanner"
